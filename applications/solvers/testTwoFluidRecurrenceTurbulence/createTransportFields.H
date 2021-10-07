/* --------------------------------------------------------------------------------- */
/* 								read flotation properties							 */
/* --------------------------------------------------------------------------------- */

Info<< "Reading scalarTransportProperties\n" << endl;
IOdictionary scalarTransportProperties
(
    IOobject
    (
        "scalarTransportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


const scalar Sc(scalarTransportProperties.lookupOrDefault<scalar>("Sc",scalar(1.0)));


const word carrierPhaseName(scalarTransportProperties.lookup("carrierPhase"));

if (carrierPhaseName != phase1.name() && carrierPhaseName != phase2.name())
{
	FatalError << "No valid carrier phase specified" << nl
		<< "Valid phase names are: " << nl
		<< phase1.name() << ",  " << phase2.name()
		<< abort(FatalError);
}

phaseModel& carrierPhase = (carrierPhaseName == phase1.name()) ? phase1 : phase2;

const word dispersePhaseName = (carrierPhaseName == phase1.name()) ? phase2.name() : phase1.name();


volScalarField& rhoCarrier = carrierPhase.thermo().rho();
volScalarField& alphaCarrier = carrierPhase;
surfaceScalarField& alphaRhoPhiCarrier = carrierPhase.alphaRhoPhi();


volScalarField contErrCarrier
(
    "contErrCarrier",
    fvc::ddt(alphaCarrier, rhoCarrier)
);



Info<< "Reading field C\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);




