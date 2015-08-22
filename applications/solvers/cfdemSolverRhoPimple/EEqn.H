// total energy equation according to
// Crowe et al.: "Multiphase flows with droplets and particles", CRC Press 1998
{
    volScalarField& he = thermo.he();
    volScalarField& T = thermo.T();
    volScalarField keff = voidfraction*turbulence->kappaEff()
            + (1-voidfraction)*particleCloud.heatTransferM().thermCondPart();
    particleCloud.heatTransferM().partFluidHeatTrans(Qsource);
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
      // + particle-fluid energy transfer due to work:
            //- fvc::div(fvc::interpolate(Us),(1-voidfraction)*p) +
      // + fluid energy dissipation due to shearing
     // net heat transfer from particles to fluid
      - Qsource
     // thermal conduction of the fluid with effective conductivity
      - fvc::laplacian(keff, T)
     ==
        fvOptions(rho, he, T)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
}
