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

// For implict T terms in the energy/enthalpy transport equation, use
// (he_n+1 - he_n) / (T_n+1 - T_n) = Cpv to eliminate T_n+1 with he_n+1.
// This formula is valid for ideal gases with e=e(T) and h=h(T). For
// incompressible fluids, e=e(T) holds, too, but enthalpy would need correction
// terms accounting for pressure variations.

    fvScalarMatrix EEqn
    (
        fvm::ddt(rhoeps, he) + fvm::div(phi, he)
      + addSource
      - Qsource
      - QCoeff*T
      - fvm::Sp(QCoeff/Cpv, he)
      + QCoeff/Cpv*he
      - fvc::laplacian(voidfraction*thCond,T)
      - fvm::laplacian(voidfraction*thCond/Cpv,he)
      + fvc::laplacian(voidfraction*thCond/Cpv,he)
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
