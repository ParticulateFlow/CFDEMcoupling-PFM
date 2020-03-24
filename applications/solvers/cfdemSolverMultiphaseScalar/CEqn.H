        // get mixture properties
		Cs  = mixture.Cs();
		diffusionCorrection = mixture.diffusionCorrection();
        Deff=particleCloud.diffCoeffM().diffCoeff(); 

        // get scalar source from DEM        
        //particleCloud.energyContributions(Qsource);
        //particleCloud.energyCoefficients(QCoeff);

        fvScalarMatrix CEqn
        (
	 	fvm::ddt(voidfraction,C)
		+ fvm::div(phi,C)
	 	- fvm::laplacian(Deff*voidfraction,C)
		+ fvm::div(fvc::interpolate(Deff*voidfraction)*diffusionCorrection*mesh.magSf(), C)
         ==
	 Csource + fvm::Sp(CCoeff,C)
        );

        CEqn.relax();
        fvOptions.constrain(CEqn);
        CEqn.solve();

        particleCloud.clockM().start(31,"postFlow");
        particleCloud.postFlow();
        particleCloud.clockM().stop("postFlow");
