fvVectorMatrix UEqn
(
    fvm::ddt(voidfractionNext,U) + MRF.DDt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
 +  (lambda*(1-voidfractionNext)/U.mesh().time().deltaT())*(fvc::Sp(1,Us)-fvm::Sp(1,U))
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
    fvOptions.correct(U);
}
