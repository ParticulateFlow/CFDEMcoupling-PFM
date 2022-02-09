volScalarField rAU = 1.0/UEqn.A();

surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));

volScalarField rAUvoidfraction("(voidfraction2|A(U))",rAU*voidfraction);

if (modelType=="A")
    rAUvoidfraction *= voidfraction;

volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(voidfraction)*fvc::flux(HbyA) );

volVectorField Uvoidfraction("(Uvoidfraction)", U * voidfraction);

surfaceScalarField phiS(fvc::interpolate(Us*voidfraction) & mesh.Sf());

phi = phiHbyA + rAUf*(fvc::interpolate(Ksl/rho) * phiS);

// rotating frames of references not tested yet
MRF.makeRelative(phi);

// adjustment of phi (only in cases w.o. p boundary conditions) not tested yet
adjustPhi(phi, U, p);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, Uvoidfraction, phiHbyA, rAUvoidfraction, MRF);


// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUvoidfraction, p) == fvc::div(phi) + particleCloud.ddtVoidfraction()
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

    if (piso.finalNonOrthogonalIter())
    {
        phi -= pEqn.flux();
    }
}

#include "continuityErrorPhiPU.H"

if (modelType=="B" || modelType=="Bfull")
    U = HbyA - rAU*fvc::grad(p) + Ksl/rho*Us*rAU;
else
    U = HbyA - voidfraction*rAU*fvc::grad(p) + Ksl/rho*Us*rAU;

U.correctBoundaryConditions();
fvOptions.correct(U);
