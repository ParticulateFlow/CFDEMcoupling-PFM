// get scalar source from DEM        
        particleCloud.forceM(1).manipulateScalarField(ETsource);
        ETsource.correctBoundaryConditions();

       // solve scalar transport equation
        fvScalarMatrix TEqn
        (
           fvm::ddt(voidfraction,T) - fvm::Sp(fvc::ddt(voidfraction),T)
         + fvm::div(phi, T) - fvm::Sp(fvc::div(phi),T)
	   - fvc::div(kT*voidfraction/(rho*Cp)*fvc::grad(T))
         ==
           ETsource/(rho*Cp)
        );

       // solve scalar transport equation
/*        fvScalarMatrix TEqn
        (
	   fvm::ddt(voidfraction*rho*Cp,T)
           + fvm::div(phi*fvc::interpolate(rho*Cp),T)
	   - fvc::div(kT*voidfraction*fvc::grad(T))
         ==
           ETsource
        );
*/


        TEqn.relax();
        TEqn.solve();
