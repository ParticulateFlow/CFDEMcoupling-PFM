// get scalar source from DEM        
particleCloud.energyContributions(ETsource);
//particleCloud.energyCoefficients(ETcoeff);


       // solve scalar transport equation
        fvScalarMatrix TEqn
        (
	 fvm::ddt(voidfraction,T) - fvm::Sp(fvc::ddt(voidfraction),T)
         + fvm::div(phi, T) - fvm::Sp(fvc::div(phi),T)
	   - fvc::div(kf*voidfraction/(rho*Cp)*fvc::grad(T))
         ==
           ETsource/(rho*Cp)
        );

/*
        fvScalarMatrix TEqn
        (
	 fvm::ddt(voidfraction*rho*Cp,T)
	 + fvm::div(phi*fvc::interpolate(rho*Cp),T)
	 - fvc::div(kf*voidfraction*fvc::grad(T))
         ==
	 ETsource //+ fvm::Sp(ETcoeff,T)
        );
*/

        TEqn.relax();
        fvOptions.constrain(TEqn);
        TEqn.solve();

        particleCloud.clockM().start(31,"postFlow");
        particleCloud.postFlow();
        particleCloud.clockM().stop("postFlow");
