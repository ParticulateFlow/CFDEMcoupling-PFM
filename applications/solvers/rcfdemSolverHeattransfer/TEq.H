	volScalarField rhoeps = rhoRec*voidfractionRec;
	
	particleCloud.energyContributions(Qsource);

        thDiff=particleCloud.thermCondM().thermDiff();
	
	dimensionedScalar Cp = particleCloud.thermCondM().Cp();
	
	fvScalarMatrix TEqn = 
	(
		fvm::ddt(rhoeps, T)
		+ fvm::div(phiRec, T)
		- fvm::laplacian(rhoeps*thDiff, T)
		- Qsource/Cp
		==
		fvOptions(rhoeps, T)	// no fvOptions support yet
	);
	
	//TEqn.relax(relaxCoeff);
	
	fvOptions.constrain(TEqn);	// no fvOptions support yet
	
	TEqn.solve();
