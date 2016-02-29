// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(voidfraction,rho, U)
   + fvm::div(phi, U)
   + particleCloud.divVoidfractionTau(U, voidfraction)
   + fvm::Sp(Ksl,U)
 ==
    fvOptions(rho, U)
);

UEqn().relax();

fvOptions.constrain(UEqn());
if (pimple.momentumPredictor() && (modelType=="B" || modelType=="Bfull"))
{
    solve(UEqn() == -fvc::grad(p)+ Ksl*Us);

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
else if (pimple.momentumPredictor())
{
    solve(UEqn() == -voidfraction*fvc::grad(p)+ Ksl*Us);

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}


