// build equation system
/*
    Note the use of the effective viscosity, which is provided by the turbulence model
    The recurrence-based turbulence models are derived from the standard base classes 
    of OpenFOAM, thus they behave as a normal turbulence model would.
*/

alphaRhoPhiCarrier = linearInterpolate(alpha2*rhoCarrier)*phi2;

fvScalarMatrix CEqn
(
	fvm::ddt(alphaCarrier*rhoCarrier, C)
  + fvm::div(alphaRhoPhiCarrier, C, "div(alphaRhoPhi,C)")
  - fvm::Sp(fvc::div(alphaRhoPhiCarrier), C)

  - fvm::laplacian
	(
		fvc::interpolate(alpha2)
	   *fvc::interpolate(carrierPhase.turbulence().muEff()/Sc),
		C
	)
	==
	fvm::SuSp(alphaCarrier*(1.0 - alphaCarrier)*rhoCarrier*K, C)
	+ fvOptions(alphaCarrier*rhoCarrier, C)
);

// solve equations
fvOptions.constrain(CEqn);
CEqn.solve();
