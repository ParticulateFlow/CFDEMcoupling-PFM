// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(voidfraction,U)
  + fvm::div(phi, U)
  + particleCloud.divVoidfractionTau(U, voidfraction)
  == 
  - fvm::Sp(Ksl/rho,U) 
);

UEqn().relax();
sources.constrain(UEqn());

if (pimple.momentumPredictor())
{
    //solve(UEqn() == -fvc::grad(p) + sources(U));
    if (modelType=="B")
        solve(UEqn() == - fvc::grad(p) + Ksl/rho*Us);
    else
        solve(UEqn() == - voidfraction*fvc::grad(p) + Ksl/rho*Us);
}
