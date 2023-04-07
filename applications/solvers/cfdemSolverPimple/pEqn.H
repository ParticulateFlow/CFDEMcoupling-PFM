volScalarField rAUvoidfraction("(voidfraction2|A(U))",rAU*voidfraction);

surfaceScalarField rAUfvoidfraction("(voidfraction2|A(U)F)", fvc::interpolate(rAUvoidfraction));

volVectorField HbyA("HbyA", U);

HbyA = rAU*UEqn.H();

phi = voidfractionf*phiByVoidfraction;

surfaceScalarField phiHbyA
(
	"phiHbyA",
	(
	    (fvc::interpolate(HbyA) & mesh.Sf() )
	    + phicForces                               //explicit contribution
	    + rAUfvoidfraction*fvc::ddtCorr(U, phiByVoidfraction) //correction
	)
);

if (modelType=="A")
    rAUvoidfraction = volScalarField("(voidfraction2|A(U))",rAU*voidfraction*voidfraction);
    
// Update the fixedFluxPressure BCs to ensure flux consistency
if (modelType=="A")
{
	volScalarField rUsed = rAU*voidfraction;
	constrainPressure(p, U, phiHbyA, rUsed,MRF);
}
else constrainPressure(p, U, phiHbyA, rAU,MRF);

while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUvoidfraction, p) == fvc::div(voidfractionf*phiHbyA) + particleCloud.ddtVoidfraction()
    );
    pEqn.setReference(pRefCell, pRefValue);
    
    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    
    if (pimple.finalNonOrthogonalIter())
    {
	    phiByVoidfraction = phiHbyA - pEqn.flux()/voidfractionf;
	    phi = voidfractionf*phiByVoidfraction;

	    #include "continuityErrorPhiPU.H"

	    // Explicitly relax pressure for momentum corrector
	    p.relax();

	    U = fvc::reconstruct(phiHbyA)
	      - rAU*fvc::reconstruct(pEqn.flux()/voidfractionf/rAUf);

	    U.correctBoundaryConditions();
	    fvOptions.correct(U);
    }
}
