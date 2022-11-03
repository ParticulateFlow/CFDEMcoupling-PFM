// Solve the Momentum equation
particleCloud.otherForces(fOther);

fvVectorMatrix UEqn
(
    fvm::div(phi, U)
  + particleCloud.divVoidfractionTau(U, voidfractionRec)
  + fvm::Sp(Ksl,U)
  - fOther
 ==
    fvOptions(rho, U)
);

if (totalStepCounter%nEveryFlow==0)
{
    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (modelType=="B" || modelType=="Bfull")
    {
        solve(UEqn == -fvc::grad(p)+ Ksl*UsRec);
    }
    else
    {
        solve(UEqn == -voidfractionRec*fvc::grad(p)+ Ksl*UsRec);
    }


    #include "limitU.H"

    fvOptions.correct(U);

    K = 0.5*magSqr(U);
}
