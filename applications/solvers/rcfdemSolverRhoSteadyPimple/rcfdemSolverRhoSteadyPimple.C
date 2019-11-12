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
    rcfdemSolverRhoSteadyPimple

Description
    Transient (DEM) + steady-state (CFD) solver for compressible flow using the
    flexible PIMPLE (PISO-SIMPLE) algorithm. Particle-motion is obtained from
    a recurrence process.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver rhoPimpleFoam in OpenFOAM(R) 4.x,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "bound.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "cfdemCloudRec.H"
#include "recBase.H"
#include "recModel.H"
#include "recPath.H"

#include "cfdemCloudEnergy.H"
#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"
#include "thermCondModel.H"
#include "energyModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createFvOptions.H"

    // create cfdemCloud
    //#include "readGravitationalAcceleration.H"
    cfdemCloudRec<cfdemCloudEnergy> particleCloud(mesh);
    #include "checkModelType.H"
    recBase recurrenceBase(mesh);
    #include "updateFields.H"

    turbulence->validate();
    //#include "compressibleCourantNo.H"
    //#include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    label recTimeIndex = 0;
    scalar recTimeStep = recurrenceBase.recM().recTimeStep();
    scalar startTime = runTime.startTime().value();

    const IOdictionary& couplingProps = particleCloud.couplingProperties();
    label nEveryFlow(couplingProps.lookupOrDefault<label>("nEveryFlow",1));
    Info << "Solving flow equations every " << nEveryFlow << " steps.\n" << endl;
    label stepcounter = 0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        particleCloud.clockM().start(1,"Global");

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // do particle stuff
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);

        //voidfraction = voidfractionRec;
        //Us = UsRec;

        if(hasEvolved)
        {
            particleCloud.smoothingM().smoothen(particleCloud.forceM(0).impParticleForces());
        }

        Info << "update Ksl.internalField()" << endl;
        Ksl = particleCloud.momCoupleM(0).impMomSource();
        Ksl.correctBoundaryConditions();

        //Force Checks
        vector fTotal(0,0,0);
        vector fImpTotal = sum(mesh.V()*Ksl.primitiveFieldRef()*(Us.primitiveFieldRef()-U.primitiveFieldRef()));
        reduce(fImpTotal, sumOp<vector>());
        Info << "TotalForceExp: " << fTotal << endl;
        Info << "TotalForceImp: " << fImpTotal << endl;

        #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");
        volScalarField rhoeps("rhoeps",rho*voidfractionRec);
        if (stepcounter%nEveryFlow==0)
        {
            while (pimple.loop())
            {
                // if needed, perform drag update here
                if (pimple.nCorrPIMPLE() <= 1)
                {
                    #include "rhoEqn.H"
                }

                // --- Pressure-velocity PIMPLE corrector loop

                #include "UEqn.H"
                #include "EEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    // besides this pEqn, OF offers a "pimple consistent"-option
                    #include "pEqn.H"
                    rhoeps=rho*voidfractionRec;
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }
        }
        stepcounter++;
        particleCloud.clockM().stop("Flow");

        particleCloud.clockM().start(31,"postFlow");
        particleCloud.postFlow();
        particleCloud.clockM().stop("postFlow");

        particleCloud.clockM().start(32,"ReadFields");
        if ( runTime.timeOutputValue() - startTime - (recTimeIndex+1)*recTimeStep + 1.0e-5 > 0.0 )
        {
            recurrenceBase.updateRecFields();
            #include "updateFields.H"
            recTimeIndex++;
        }
        particleCloud.clockM().stop("ReadFields");

        runTime.write();


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;


        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
