fvVectorMatrix UEqn
(
    fvm::ddt(voidfraction,U) + MRF.DDt(U)
  // + 1.5*fvm::div(phi, U)                //Adams Bashforth implicit scheme
  // - 0.5*fvc::div(phi, U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
  + V_ratio*lambda*Us/U.mesh().time().deltaT()
  - fvm::Sp(V_ratio*lambda/U.mesh().time().deltaT(),U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
    fvOptions.correct(U);
}
