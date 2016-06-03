        // get scalar source from DEM        
        particleCloud.energyContributions(Qsource);
	
        // solve temperature transport equation assuming const. density and heat capacity
	
 	thDiff=particleCloud.thermCondM().thermDiff();
        // BC for thDiff are calculated from its mathematical expression
        // to enforce user-defined values, change thDiff to MUST_READ and do
        // thDiff.correctBoundaryConditions();
		
        fvScalarMatrix TEqn
        (
           fvm::ddt(voidfraction,T) + fvm::div(phi, T)
         - fvm::Sp(fvc::ddt(voidfraction),T) - fvm::Sp(fvc::div(phi),T)
	 - fvm::laplacian(voidfraction*thDiff, T)
         ==
           Qsource/(rho*Cp)
        );
	
	    // work in progress:
    // if (particleCloud.implicitEnergyModel)
    // volScalarField& tpart (lookup ...)
    // particleCloud.energyCoefficients(Qcoeff);
    // TEqn += ...
	
       // TEqn.relax();
        TEqn.solve();
	
	gradT=fvc::grad(T);