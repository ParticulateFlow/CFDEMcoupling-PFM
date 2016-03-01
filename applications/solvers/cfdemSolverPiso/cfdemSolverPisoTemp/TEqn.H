        // get scalar source from DEM        
        particleCloud.energyContributions(Qsource);
	
        // solve temperature transport equation assuming const. density and heat capacity
        fvScalarMatrix TEqn
        (
           fvm::ddt(voidfraction,T) - fvm::Sp(fvc::ddt(voidfraction),T)
         + fvm::div(phi, T) - fvm::Sp(fvc::div(phi),T)
         - fvm::laplacian(voidfraction*particleCloud.thermCondM().thermDiff(), T)
         ==
           Qsource/(rho*Cp)
        );
        TEqn.relax();
        TEqn.solve();
