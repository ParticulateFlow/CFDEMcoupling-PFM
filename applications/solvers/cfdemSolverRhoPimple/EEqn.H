// contributions to internal energy equation can be found in
// Crowe et al.: "Multiphase flows with droplets and particles", CRC Press 1998
//
// fluid thermal conductivity according to
// Syamlal, Gidaspow. AIChE J. 31 (1985)
{
    volScalarField& he = thermo.he();
    volScalarField keff = (1-Foam::sqrt(1-voidfraction+SMALL)) / (voidfraction+SMALL);
    keff *= turbulence->kappaEff();
    particleCloud.energyM(0).energyContribution(Qsource);
    Qsource.correctBoundaryConditions();

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
      - fvc::laplacian(keff, T)
      // + particle-fluid energy transfer due to work
      // + fluid energy dissipation due to shearing
     ==
        fvOptions(rho, he, T)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
}
