        // get scalar source from DEM        
        particleCloud.energyContributions(Qsource);
	
        // solve temperature transport equation assuming const. density and heat capacity
	
	thDiff=voidfraction*particleCloud.thermCondM().thermDiff();
	thDiff.correctBoundaryConditions();
		
        fvScalarMatrix TEqn
        (
           fvm::ddt(voidfraction,T) + fvm::div(phi, T)
         - fvm::Sp(fvc::ddt(voidfraction),T) - fvm::Sp(fvc::div(phi),T)
       //  - fvm::laplacian(voidfraction*particleCloud.thermCondM().thermDiff(), T)
	 - fvm::laplacian(thDiff, T)
         ==
           Qsource/(rho*Cp)
        );
	
       // TEqn.relax();
        TEqn.solve();
	
	gradT=fvc::grad(T);