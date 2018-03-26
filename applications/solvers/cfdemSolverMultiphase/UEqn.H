const surfaceScalarField& rhoPhi(mixture.rhoPhi());

fvVectorMatrix UEqn
(
    fvm::ddt(rhoEps, U) - fvm::Sp(fvc::ddt(rhoEps),U)
  + fvm::div(rhoPhi, U) - fvm::Sp(fvc::div(rhoPhi),U)
  + particleCloud.divVoidfractionTau(U, rhoEps)
 ==
    fvOptions(rho, U)
  - fvm::Sp(Ksl,U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor() && (modelType=="B" || modelType=="Bfull"))
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (- ghf*fvc::snGrad(rho) - fvc::snGrad(p_rgh)) * mesh.magSf()
        )
     +
        fvc::reconstruct
        (
            mixture.surfaceTensionForce() * mesh.magSf()
        ) * voidfraction
     + Ksl*Us
    );
    
    fvOptions.correct(U);
}        
else if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                mixture.surfaceTensionForce()
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        ) * voidfraction
     + Ksl*Us
    );

    fvOptions.correct(U);
}
