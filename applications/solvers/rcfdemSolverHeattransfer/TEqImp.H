	volScalarField rhoeps = rhoRec*voidfractionRec;
	
	particleCloud.energyContributions(Qsource);

        thDiff=particleCloud.thermCondM().thermDiff();
	
	dimensionedScalar Cp = particleCloud.thermCondM().Cp()*718/1005; // need cV, hence the factor

	K = 0.5*magSqr(URec);

	addSource = fvc::ddt(rhoeps, K) + fvc::div(phiRec, K) +
	fvc::div
	(
            fvc::absolute(phiRec/fvc::interpolate(rhoRec), voidfractionRec*URec), pRec
	);

	particleCloud.energyCoefficients(QCoeff);
	
	fvScalarMatrix TEqn = 
	(
		fvm::ddt(rhoeps, T)
		+ fvm::div(phiRec, T)
		+ addSource/Cp
		- fvm::laplacian(rhoeps*thDiff, T)
		- Qsource/Cp
		- fvm::Sp(QCoeff/Cp, T)
		==
		fvOptions(rhoeps, T)	// no fvOptions support yet
	);
	
	//TEqn.relax(relaxCoeff);
	
	fvOptions.constrain(TEqn);	// no fvOptions support yet
	
	TEqn.solve();

        particleCloud.clockM().start(31,"postFlow");
	counter++;
	if( (counter - couplingSubStep) % dtDEM2dtCFD == 0)
	    particleCloud.postFlow();
        particleCloud.clockM().stop("postFlow");

	
	// change couplingSubStep during course of simulation
//	if( counter == 16000 )
//	    couplingSubStep++;
	
//	if (doCouple)
//	    particleCloud.postFlow();
