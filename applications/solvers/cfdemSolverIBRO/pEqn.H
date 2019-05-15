volScalarField rUA = 1.0/UEqn.A();
surfaceScalarField rUAf(fvc::interpolate(rUA));

U = rUA*UEqn.H();

surfaceScalarField phiS(fvc::interpolate(Us*lambda) & mesh.Sf());

phi = (fvc::interpolate(U) & mesh.Sf())
    + rUAf*fvc::ddtCorr(U, phi);            // Is there additional flux term due to the particle presence?

MRF.makeRelative(phi);

adjustPhi(phi, U, p);

while (piso.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rUA, p) == fvc::div(phi) + particleCloud.ddtVoidfraction()
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

    if (piso.finalNonOrthogonalIter())
    {
        phi -= pEqn.flux();
    }
}

#include "continuityErrs.H"
U -= rUA*fvc::grad(p); // should we add a pressure correction?
U.correctBoundaryConditions();
fvOptions.correct(U);
