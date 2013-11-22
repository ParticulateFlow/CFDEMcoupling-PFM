// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(voidfraction,U)
  + fvm::div(phi, U)
  + particleCloud.divVoidfractionTau(U, voidfraction)
  == 
  - fvm::Sp(Ksl/rho,U)
  #if defined(version22)
  + fvOptions(U)  // test phase
  #endif
);

UEqn().relax();
#if defined(version22)
fvOptions.constrain(UEqn());
#endif
//sources.constrain(UEqn());

if (pimple.momentumPredictor())
{
    //solve(UEqn() == -fvc::grad(p) + sources(U));
    if (modelType=="B")
        solve(UEqn() == - fvc::grad(p) + Ksl/rho*Us);
    else
        solve(UEqn() == - voidfraction*fvc::grad(p) + Ksl/rho*Us);

    #if defined(version22)
    fvOptions.correct(U);  // test phase
    #endif
}
