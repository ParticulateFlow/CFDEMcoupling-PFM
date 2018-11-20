rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rhoeps*rAU));
if (modelType=="A")
{
    rhorAUf *= fvc::interpolate(voidfractionRec);
}
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
//tUEqn.clear();

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

    while (pimple.correctNonOrthogonal())
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvc::div(phi)
          - fvm::laplacian(rhorAUf, p)
          ==
            fvOptions(psi, p, rho.name())
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pEqn.flux();
        }
    }
}

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
    U = HbyA - rAU*(voidfractionRec*fvc::grad(p)-Ksl*UsRec);
}
else
{
    U = HbyA - rAU*(fvc::grad(p)-Ksl*UsRec);
}

U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);
