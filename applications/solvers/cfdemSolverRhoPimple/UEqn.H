// Solve the Momentum equation
particleCloud.otherForces(fOther);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rhoeps, U)
  + fvm::div(phi, U)
  + particleCloud.divVoidfractionTau(U, voidfraction)
  + fvm::Sp(Ksl,U)
  - fOther
 ==
    fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor() && (modelType=="B" || modelType=="Bfull"))
{
    solve(UEqn == -fvc::grad(p)+ Ksl*Us);

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
else if (pimple.momentumPredictor())
{
    solve(UEqn == -voidfraction*fvc::grad(p)+ Ksl*Us);

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
