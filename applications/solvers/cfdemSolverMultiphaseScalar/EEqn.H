        // get mixture properties
        Cp = mixture.Cp();
        kf = mixture.kf();

        // get scalar source from DEM
        particleCloud.energyContributions(Qsource);
        particleCloud.energyCoefficients(QCoeff);

        fvScalarMatrix EEqn
        (
            rho*Cp*(fvm::ddt(voidfraction,T)
          + fvm::div(phi,T))
          - fvm::laplacian(thCond*voidfraction,T)
          ==
            Qsource + fvm::Sp(QCoeff,T)
        );


        EEqn.relax();
        fvOptions.constrain(EEqn);
        EEqn.solve();

