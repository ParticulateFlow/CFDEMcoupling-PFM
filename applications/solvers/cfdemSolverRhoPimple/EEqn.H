// contributions to internal energy equation can be found in
// Crowe et al.: "Multiphase flows with droplets and particles", CRC Press 1998
{
    // dim he = J / kg
    volScalarField& he = thermo.he();
    particleCloud.energyContributions(Qsource);

    thDiff=particleCloud.thermCondM().thermDiff();   
    // BC for thDiff are calculated from its mathematical expression
    // to enforce user-defined values, change thDiff to MUST_READ and do
    // thDiff.correctBoundaryConditions();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rhoeps, he) + fvm::div(phi, he)
      + fvc::ddt(rhoeps, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), voidfraction*U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
     // net heat transfer from particles to fluid
      - Qsource
     // thermal conduction of the fluid with effective conductivity
      - fvm::laplacian(rhoeps*thDiff,he)
          
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
}
