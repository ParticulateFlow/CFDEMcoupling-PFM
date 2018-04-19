// Solve the Momentum equation
particleCloud.otherForces(fOther);

tmp<fvVectorMatrix> tUEqn
(
    fvm::div(phi, U)
  + particleCloud.divVoidfractionTau(U, voidfraction)
  + fvm::Sp(Ksl,U)
  - fOther
 ==
    fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (modelType=="B" || modelType=="Bfull")
{
    solve(UEqn == -fvc::grad(p)+ Ksl*Us);
}
else
{
    solve(UEqn == -voidfraction*fvc::grad(p)+ Ksl*Us);
}

fvOptions.correct(U);

// limit vel

forAll(U,cellI)
{
    scalar mU = mag(U[cellI]);
    if (mU > 200.0) U[cellI] *= 200.0/mU;
}


K = 0.5*magSqr(U);
