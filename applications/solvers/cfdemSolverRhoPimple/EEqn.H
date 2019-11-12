// contributions to internal energy equation can be found in
// Crowe et al.: "Multiphase flows with droplets and particles", CRC Press 1998
{
    // dim he = J / kg
    volScalarField& he = thermo.he();
    particleCloud.energyContributions(Qsource);
    particleCloud.energyCoefficients(QCoeff);

    addSource = fvc::ddt(rhoeps, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), voidfraction*U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        );

    Cpv = he.name() == "e" ? thermo.Cv() : thermo.Cp();

    // correct source for the thermodynamic reference temperature
    dimensionedScalar Tref("Tref", dimTemperature, T[0]-he[0]/(Cpv[0]+SMALL));
    Qsource += QCoeff*Tref;

    fvScalarMatrix EEqn
    (
        fvm::ddt(rhoeps, he) + fvm::div(phi, he)
      + addSource
     // net heat transfer from particles to fluid
      - Qsource
      - fvm::Sp(QCoeff/Cpv, he)
     // thermal conduction of the fluid with effective conductivity
      - fvm::laplacian(voidfraction*thCond/Cpv,he)
      // + particle-fluid energy transfer due to work
      // + fluid energy dissipation due to shearing
     ==
        fvOptions(rho, he)
    );


    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "T max/min : " << max(T).value() << " " << min(T).value() << endl;

    particleCloud.clockM().start(31,"energySolve");
    particleCloud.solve();
    particleCloud.clockM().stop("energySolve");
}
