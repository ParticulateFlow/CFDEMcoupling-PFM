/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2015- Thomas Lichtenegger, JKU Linz, Austria

Application
    cfdemSolverRhoPimpleChem

Description
    Transient solver for compressible flow using the flexible PIMPLE (PISO-SIMPLE)
    algorithm.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver rhoPimpleFoam in OpenFOAM(R) 2.3,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulenceModel.H"
#include "psiCombustionModel.H"
#include "bound.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"

#include "cfdemCloudEnergy.H"
#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"
#include "thermCondModel.H"
#include "energyModel.H"
#include "chemistryModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    
    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    cfdemCloudEnergy particleCloud(mesh);
    #include "checkModelType.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
	
	particleCloud.clockM().start(1,"Global");

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // If DEM particle does not exist (because radius_[i] < rmin, stop coupling!
        // do particle stuff
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);

        if(hasEvolved)
        {
            particleCloud.smoothingM().smoothen(particleCloud.forceM(0).impParticleForces());
        }
    
        Info << "update Ksl.internalField()" << endl;
        Ksl = particleCloud.momCoupleM(0).impMomSource();
        Ksl.correctBoundaryConditions();

       //Force Checks
       vector fTotal(0,0,0);
       vector fImpTotal = sum(mesh.V()*Ksl.internalField()*(Us.internalField()-U.internalField()));
       reduce(fImpTotal, sumOp<vector>());
       Info << "TotalForceExp: " << fTotal << endl;
       Info << "TotalForceImp: " << fImpTotal << endl;

       Info << "number of particles = " << particleCloud.numberOfParticles() << nl << endl;

        #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");



        particleCloud.clockM().start(26,"Flow");
	
        if (pimple.nCorrPIMPLE() <= 1)
        {
            #include "rhoEqn.H"
        }
        volScalarField rhoeps("rhoeps",rho*voidfraction);
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "YEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
	      	rhoeps=rho*voidfraction;
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

       
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
