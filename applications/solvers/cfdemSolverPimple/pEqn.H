volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);

surfaceScalarField rUAfvoidfraction("(voidfraction2|A(U)F)", fvc::interpolate(rUAvoidfraction));

volVectorField HbyA("HbyA", U);

HbyA = rUA*UEqn.H();

phi = voidfractionf*phiByVoidfraction;

surfaceScalarField phiHbyA
(
	"phiHbyA",
	(
	    (fvc::interpolate(HbyA) & mesh.Sf() )
	    + phicForces                               //explicit contribution
	    + rUAfvoidfraction*fvc::ddtCorr(U, phiByVoidfraction) //correction
	)
);

if (modelType=="A")
    rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);
    
// Update the fixedFluxPressure BCs to ensure flux consistency
if (modelType=="A")
{
	volScalarField rUsed = rUA*voidfraction;
	constrainPressure(p, U, phiHbyA, rUsed,MRF);
}
else constrainPressure(p, U, phiHbyA, rUA,MRF);

while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rUAvoidfraction, p) == fvc::div(voidfractionf*phiHbyA) + particleCloud.ddtVoidfraction()
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
	      - rUA*fvc::reconstruct(pEqn.flux()/voidfractionf/rUAf);

	    U.correctBoundaryConditions();
	    fvOptions.correct(U);
    }
}
