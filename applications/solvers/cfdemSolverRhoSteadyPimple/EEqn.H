// contributions to internal energy equation can be found in
// Crowe et al.: "Multiphase flows with droplets and particles", CRC Press 1998
{
    // dim he = J / kg
    volScalarField& he = thermo.he();
    particleCloud.energyContributions(Qsource);
    particleCloud.energyCoefficients(QCoeff);

    //thDiff=particleCloud.thermCondM().thermDiff();   
    thCond=particleCloud.thermCondM().thermCond(); 

    addSource = 
        (
            he.name() == "e"
          ? 
            fvc::div(phi, K) +
            fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), voidfraction*U),
                p,
                "div(phiv,p)"
            )
          : fvc::div(phi, K)
        );

    Cpv = he.name() == "e" ? thermo.Cv() : thermo.Cp();


    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + addSource
      - Qsource
      - fvm::Sp(QCoeff/Cpv, he)
      - fvm::laplacian(voidfraction*thCond/Cpv,he)
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
