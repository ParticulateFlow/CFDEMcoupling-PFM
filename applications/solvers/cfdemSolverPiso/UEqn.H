fvVectorMatrix UEqn
(
    fvm::ddt(voidfraction,U) - fvm::Sp(fvc::ddt(voidfraction),U)
  + fvm::div(phi,U) - fvm::Sp(fvc::div(phi),U)
  + particleCloud.divVoidfractionTau(U, voidfraction)
  ==
    fvOptions(U)
  - fvm::Sp(Ksl/rho,U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (piso.momentumPredictor() && (modelType=="B" || modelType=="Bfull"))
{
    solve(UEqn == - fvc::grad(p) + Ksl/rho*Us);
    fvOptions.correct(U);
}
else if (piso.momentumPredictor())
{
    solve(UEqn == - voidfraction*fvc::grad(p) + Ksl/rho*Us);
    fvOptions.correct(U);
}
