        // get mixture properties
        D  = mixture.D();
		Cs  = mixture.Cs();
		diffusionCorrection = mixture.diffusionCorrection();
        //keff=particleCloud.thermCondM().thermCond(); 

        // get scalar source from DEM        
        //particleCloud.energyContributions(Qsource);
        //particleCloud.energyCoefficients(QCoeff);

        fvScalarMatrix CEqn
        (
	 	fvm::ddt(voidfraction,C)
		+ fvm::div(phi,C)
	 	- fvm::laplacian(D*voidfraction,C)
		//+ fvm::div(fvc::interpolate(D*voidfraction*diffusionCorrection) & mesh.Sf(), C)
         ==
	 Csource + fvm::Sp(CCoeff,C)
        );

        CEqn.relax();
        fvOptions.constrain(CEqn);
        CEqn.solve();

		//diffusiveFlux = fvc::laplacian(D*voidfraction,C);
		//correctionFlux = fvc::div(D*voidfraction*diffusionCorrection*C);

        particleCloud.clockM().start(31,"postFlow");
        particleCloud.postFlow();
        particleCloud.clockM().stop("postFlow");
