{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(voidfraction,rho)
      + fvc::div(phi)
      ==
        particleCloud.chemistryM(0).Sm()
      + fvOptions(rho)
    );

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho);
}

// ************************************************************************* //
