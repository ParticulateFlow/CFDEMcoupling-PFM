/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverPisoMS

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6, 
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"

#include "cfdemCloudMS.H"
#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    cfdemCloudMS particleCloud(mesh);
    #include "checkModelType.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (runTime.loop())
    {
        Info<< "\nStarting time loop\n" << endl;
            particleCloud.clockM().start(1,"Global");

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // do particle stuff
        particleCloud.clockM().start(2,"Coupling");
        particleCloud.evolve(voidfraction,Us,U);
      
        Info << "update Ksl.internalField()" << endl;
        Ksl = particleCloud.momCoupleM(0).impMomSource();
        particleCloud.smoothingM().smoothen(Ksl);
        Ksl.correctBoundaryConditions();

        #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");
        // Pressure-velocity PISO corrector
        {
            // Momentum predictor
            fvVectorMatrix UEqn
            (
                fvm::ddt(voidfraction,U) //particleCloud.ddtVoidfractionU(U,voidfraction) //
              + fvm::div(phi, U)
//              + turbulence->divDevReff(U)
              + particleCloud.divVoidfractionTau(U, voidfraction)
             ==
              - fvm::Sp(Ksl/rho,U)
            );

            if (modelType=="B")
                UEqn == - fvc::grad(p) + Ksl/rho*Us;
            else
                UEqn == - voidfraction*fvc::grad(p) + Ksl/rho*Us;

            UEqn.relax();

            if (piso.momentumPredictor())
                solve(UEqn);

            // --- PISO loop

            while (piso.correct())
            {
                volScalarField rUA = 1.0/UEqn.A();

                surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
                volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);

                U = rUA*UEqn.H();

                phi = (fvc::interpolate(U*voidfraction) & mesh.Sf() );
                     //+ fvc::ddtPhiCorr(rUAvoidfraction, U, phi);
                surfaceScalarField phiS(fvc::interpolate(Us*voidfraction) & mesh.Sf());
                surfaceScalarField phiGes = phi + rUAf*(fvc::interpolate(Ksl/rho) * phiS);

                if (modelType=="A")
                    rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);

                // Non-orthogonal pressure corrector loop
                while (piso.correctNonOrthogonal())
                {
                    // Pressure corrector
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rUAvoidfraction, p) == fvc::div(phiGes) + particleCloud.ddtVoidfraction()
                    );
                    pEqn.setReference(pRefCell, pRefValue);

                    pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                    if (piso.finalNonOrthogonalIter())
                    {
                        phiGes -= pEqn.flux();
                    }

                } // end non-orthogonal corrector loop

                #include "continuityErrorPhiPU.H"

                if (modelType=="B")
                    U -= rUA*fvc::grad(p) - Ksl/rho*Us*rUA;
                else
                    U -= voidfraction*rUA*fvc::grad(p) - Ksl/rho*Us*rUA;

                U.correctBoundaryConditions();

            } // end piso loop
        }

        turbulence->correct();

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
