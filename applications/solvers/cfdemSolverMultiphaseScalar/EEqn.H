        // get mixture properties
        Cp  = mixture.Cp();
		kf  = mixture.kf();
        keff=particleCloud.thermCondM().thermCond(); 

        // get scalar source from DEM        
        particleCloud.energyContributions(Qsource);
        particleCloud.energyCoefficients(QCoeff);

        fvScalarMatrix EEqn
        (
	 rho*Cp*(fvm::ddt(voidfraction,T)
		 + fvm::div(phi,T))
	 - fvc::div(keff*voidfraction*fvc::grad(T))
         ==
	 Qsource + fvm::Sp(QCoeff,T)
        );


        EEqn.relax();
        fvOptions.constrain(EEqn);
        EEqn.solve();

        particleCloud.clockM().start(31,"postFlow");
        particleCloud.postFlow();
        particleCloud.clockM().stop("postFlow");
