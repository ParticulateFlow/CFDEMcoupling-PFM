rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rhoeps*rAU));
if (modelType=="A")
{
    rhorAUf *= fvc::interpolate(voidfraction);
}
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

surfaceScalarField phiUs("phiUs", fvc::interpolate(rhoeps*rAU*Ksl*Us)& mesh.Sf());

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

if (pimple.transonic())
{
//  transonic version not implemented yet
}
else
{
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rhoeps*HbyA)
    //      + rhorAUf*fvc::ddtCorr(rho, U, phi)
        )
    );

    // flux without pressure gradient contribution
    phi = phiHbyA + phiUs;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rhoeps, U, phi, rhorAUf);

    while (pimple.correctNonOrthogonal())
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi*voidfraction, p)
          + fvc::div(phi)
          - fvm::laplacian(rhorAUf, p)
          ==
            fvOptions(psi, p, rho.name())
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrsPU.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Recalculate density from the relaxed pressure
rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
Info<< "rho max/min : " << max(rho).value()
    << " " << min(rho).value() << endl;

if (modelType=="A")
{
    U = HbyA - rAU*(voidfraction*fvc::grad(p)-Ksl*Us);
}
else
{
    U = HbyA - rAU*(fvc::grad(p)-Ksl*Us);
}
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(voidfraction,p);
}
