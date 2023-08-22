dimensionedVector gN = mesh.lookupObject<uniformDimensionedVectorField> ("g");
dimensionedScalar volume = fvc::domainIntegrate(unity);

Info<< "particle_ENSTROPHY: "
    << (
            fvc::domainIntegrate( 0.5*magSqr(fvc::curl(Us)))
            /volume
       ).value()
    << endl;

Info<< "air_ENSTROPHY: "
    << (
           fvc::domainIntegrate( 0.5*magSqr(fvc::curl(U)))
           /volume
       ).value()
    << endl;

Info<< "slip_velocity: "
    << - ((
		fvc::domainIntegrate(voidfraction*(U&gN)).value()
	       /fvc::domainIntegrate(voidfraction*mag(gN)).value()
	     )
	   - (
		fvc::domainIntegrate((1.0-voidfraction)*(Us&gN)).value()
	       /fvc::domainIntegrate((1.0-voidfraction)*mag(gN)).value()
         ))
   << endl;
   
dimensionedVector alpha1Us = fvc::domainIntegrate((1.0-voidfraction)*(Us))/fvc::domainIntegrate((1.0-voidfraction));
dimensionedVector alpha2U = fvc::domainIntegrate(voidfraction*(U))/fvc::domainIntegrate(voidfraction);
dimensionedScalar alpha1M  = fvc::domainIntegrate((1.0-voidfraction))/volume;
dimensionedScalar alpha2M  = scalar(1.0) - alpha1M;

Info<< "TKE gas: "
    << 0.5
      *(
          fvc::domainIntegrate(voidfraction*(U&U)).value()
         /fvc::domainIntegrate(voidfraction).value()
       )
     - 0.5
      *(
          alpha2U.value()
         &alpha2U.value()
       )
    << endl;

Info<< "TKE solid: "
    << 0.5
      *(
          fvc::domainIntegrate((1.0-voidfraction)*(Us&Us)).value()
         /fvc::domainIntegrate(1.0-voidfraction).value()
       )
     - 0.5
      *(
          alpha1Us.value()
         &alpha1Us.value()
       )
    << endl;

Info<< "PhiP2: "
    << fvc::domainIntegrate((1.0-voidfraction)*(1.0-voidfraction)).value()
       /fvc::domainIntegrate(unity).value()
       - alpha1M.value()*alpha1M.value()
    << endl;

