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
    cfdemSolverIB

Description
    Transient solver for incompressible flow.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling using immersed body
    (fictitious domain) method is added.
Contributions
    Alice Hager
\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"

#include "cfdemCloudIB.H"
#include "implicitCouple.H"

#include "averagingModel.H"
#include "regionModel.H"
#include "voidFractionModel.H"

#include "dynamicFvMesh.H" //dyM

#include "cellSet.H"

#include "fvOptions.H"  // fvOptions library

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    //#include "createMesh.H"
    #include "createDynamicFvMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    // create cfdemCloud for particle data
    cfdemCloudIB particleCloud(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info << "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        //=== dyM =====================
        interFace = mag(mesh.lookupObject<volScalarField>("voidfractionNext"));
        mesh.update(); //dyM

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // particle information update
        Info << "- evolve()" << endl;
        particleCloud.evolve(Us);

        // Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            // Momentum predictor

            tmp<fvVectorMatrix> tUEqn
            (
                fvm::ddt(voidfraction,U)
              + fvm::div(phi,U)
              + turbulence->divDevReff(U)
              ==
                fvOptions(U)
            );

            fvVectorMatrix UEqn = tUEqn.ref();

            UEqn.relax();

            fvOptions.constrain(UEqn);

            solve(UEqn == -fvc::grad(p));

            fvOptions.correct(U);

             // Pressure corrector loop
            while (pimple.correct())
            {
                volScalarField rAU(1.0/UEqn.A());
                volVectorField HbyA("HbyA", U);
                HbyA = rAU*UEqn.H();

                surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    fvc::flux(HbyA)
                  + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
                );

                adjustPhi(phiHbyA, U, p);

                tmp<volScalarField> rAtU(rAU);

                if (pimple.consistent())
                {
                    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
                    phiHbyA +=
                        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
                    HbyA -= (rAU - rAtU())*fvc::grad(p);
                }

                if (pimple.nCorrPISO() <= 1)
                {
                    tUEqn.clear();
                }

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(p, U, phiHbyA, rAtU());

                // Non-orthogonal pressure corrector loop
                while (pimple.correctNonOrthogonal())
                {
                    // Pressure corrector
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

                    if (pimple.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                #include "continuityErrs.H"

                p.relax();

                // Momentum corrector
                U = HbyA - rAtU()*fvc::grad(p);
                U.correctBoundaryConditions();
                fvOptions.correct(U);
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
       }

            Info << "particleCloud.calcVelocityCorrection() " << endl;
            volScalarField voidfractionNext = mesh.lookupObject<volScalarField>("voidfractionNext");
            particleCloud.calcVelocityCorrection(p,U,phiIB,voidfractionNext);

            fvOptions.correct(U);

            runTime.write();

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
   }

       Info << "End\n" << endl;

       return 0;
}

// ************************************************************************* //
