volScalarField rAU(1.0/UEqn().A());
surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rAU));

//U = rAU*(UEqn() == sources(U))().H(); // new formulation
U = rAU*UEqn().H();   // old

if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}

phi = (fvc::interpolate(U*voidfraction) & mesh.Sf())
    + fvc::ddtPhiCorr(rAU, U, phi);

//adjustPhi(phi, U, p);//??

surfaceScalarField phiS(fvc::interpolate(Us*voidfraction) & mesh.Sf());
surfaceScalarField phiGes = phi + rUAf*(fvc::interpolate(Ksl/rho) * phiS);

volScalarField rUAvoidfraction("(voidfraction2|A(U))",rAU*voidfraction);
if (modelType=="A")
    rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rAU*voidfraction*voidfraction);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rUAvoidfraction, p) == fvc::div(phiGes) //+ particleCloud.ddtVoidfraction()
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phiGes -= pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

if (modelType=="B")
    U -= rAU*fvc::grad(p) - Ksl/rho*Us*rAU;
else
    U -= voidfraction*rAU*fvc::grad(p) - Ksl/rho*Us*rAU;

U.correctBoundaryConditions();
#if defined(version22)
fvOptions.correct(U);  // test phase
#endif
//sources.correct(U);
