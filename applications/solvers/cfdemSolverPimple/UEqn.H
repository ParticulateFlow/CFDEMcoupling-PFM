particleCloud.otherForces(fOther);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(voidfraction,U) - fvm::Sp(fvc::ddt(voidfraction),U)
  + fvm::div(phi,U) - fvm::Sp(fvc::div(phi),U) 
  - fOther/rho
  ==
    fvOptions(U)
  - fvm::Sp(Ksl/rho,U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

volScalarField rUA = 1.0/UEqn.A();

surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);

surfaceScalarField phicForces
(
	fvc::interpolate(rUA*(Ksl*Us)/rho) & mesh.Sf()
);

if (pimple.momentumPredictor() && (modelType=="B" || modelType=="Bfull"))
{
    solve(UEqn == fvc::reconstruct(phicForces/rUAf - fvc::snGrad(p)*mesh.magSf()));
    fvOptions.correct(U);
}
else if (pimple.momentumPredictor())
{
    solve(UEqn == fvc::reconstruct(phicForces/rUAf - fvc::snGrad(p)*voidfractionf*mesh.magSf()));
    fvOptions.correct(U);
}

