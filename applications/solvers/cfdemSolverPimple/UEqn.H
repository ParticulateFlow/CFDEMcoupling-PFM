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

volScalarField rAU = 1.0/UEqn.A();

surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));

surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);

surfaceScalarField phicForces
(
    fvc::interpolate(rAU*(Ksl*Us)/rho) & mesh.Sf()
);

if (pimple.momentumPredictor() && (modelType=="B" || modelType=="Bfull"))
{
    solve(UEqn == fvc::reconstruct(phicForces/rAUf - fvc::snGrad(p)*mesh.magSf()));
    fvOptions.correct(U);
}
else if (pimple.momentumPredictor())
{
    solve(UEqn == fvc::reconstruct(phicForces/rAUf - fvc::snGrad(p)*voidfractionf*mesh.magSf()));
    fvOptions.correct(U);
}

