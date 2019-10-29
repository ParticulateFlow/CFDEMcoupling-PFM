fvVectorMatrix UEqn
(
    fvm::ddt(voidfractionNext,U) + MRF.DDt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
  + lambda*(1-voidfractionNext)*Us/U.mesh().time().deltaT() // explicit part of source term
  - fvm::Sp(lambda*(1-voidfractionNext)/U.mesh().time().deltaT(),U) //fvm::ddt(V_ratio*lambda,U) //   // implicit part of source term
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
    fvOptions.correct(U);
}
