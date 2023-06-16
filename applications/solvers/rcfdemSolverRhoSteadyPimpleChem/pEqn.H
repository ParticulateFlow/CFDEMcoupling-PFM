rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();

if (totalStepCounter%nEveryFlow==0)
{

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rhoeps*rAU));
if (modelType=="A")
{
    rhorAUf *= fvc::interpolate(voidfractionRec);
}
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

surfaceScalarField phiUs("phiUs", fvc::interpolate(rhoeps*rAU*Ksl*UsRec)& mesh.Sf());


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
        )
    );

    // flux without pressure gradient contribution
    phi = phiHbyA + phiUs;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rhoeps, U, phi, rhorAUf);

    volScalarField Sm0(neg(Sm)*Sm/(p + psmall));
    volScalarField Sm1(pos0(Sm)*Sm);

    while (pimple.correctNonOrthogonal())
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvc::div(phi)
          - fvm::laplacian(rhorAUf, p)
          ==
            fvm::Sp(Sm0,p)
          + Sm1
          + fvOptions(psi, p, rho.name())
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pEqn.flux();
        }
    }
    #include "monitorMassSource.H"
}

#include "rhoEqn.H"
#include "compressibleContinuityErrsPU.H"

// Explicitly relax pressure for momentum corrector
p.relax();

#include "limitP.H"

// Recalculate density from the relaxed pressure
rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
Info<< "rho max/min : " << max(rho).value()
    << " " << min(rho).value() << endl;

if (modelType=="A")
{
    U = HbyA - rAU*(voidfractionRec*fvc::grad(p)-Ksl*UsRec);
}
else
{
    U = HbyA - rAU*(fvc::grad(p)-Ksl*UsRec);
}

#include "limitU.H"

U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

}
