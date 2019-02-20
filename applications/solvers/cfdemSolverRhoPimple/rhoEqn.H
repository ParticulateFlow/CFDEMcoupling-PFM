{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(voidfraction,rho)
      + fvc::div(phi)
      ==
        fvOptions(rho)
    );

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho);
}

// ************************************************************************* //
