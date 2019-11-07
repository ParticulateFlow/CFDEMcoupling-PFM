    // dummy fields
    Info<< "\nCreating dummy pressure and density fields\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("p", dimensionSet(1, 2, -2, 0, 0), 1.0)
    );
    
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimensionSet(1, -3, 0, 0, 0), 1.0)
    );
    
    // recurrence fields
    Info<< "\nCreating recurrence fields.\n" << endl;

    volVectorField URec
    (
        IOobject
        (
            "URec",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField U2Rec
    (
        IOobject
        (
            "U2Rec",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // calculated fields
    Info<< "\nCreating fields subject to calculation\n" << endl;

    volScalarField delta
    (
        IOobject
        (
            "delta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
       mesh,
       dimensionedScalar("delta", dimLength, 0.0)
    );

    delta.primitiveFieldRef()=pow(mesh.V(),1.0/3.0);
    delta.write();


    Info<< "\ncreating dU2\n" << endl;

    volScalarField dU2
    (
        IOobject
        (
            "dU2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt(0.5*mag(U2Rec - magSqr(URec)))*delta*0.094
     );

    forAll(dU2, cellI)
    {
        if (U2Rec[cellI]-magSqr(URec[cellI]) < 0.0)
        {
            dU2[cellI] = 0.0;
        }
    }

    dU2.write();

Info<< "Calculating face flux field phiRec\n" << endl;
surfaceScalarField phiRec
 (
     IOobject
     (
        "phiRec",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
     linearInterpolate(URec) & mesh.Sf()
 );

    phiRec.write();

    singlePhaseTransportModel laminarTransport(URec, phiRec);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(URec, phiRec, laminarTransport)
    );
  
    dimensionedScalar Sc("Sc", dimless, laminarTransport);
    dimensionedScalar Sct("Sct", dimless, laminarTransport);
  
 // create concentration field
 Info<< "Creating scalar transport field\n" << endl;

    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    fvScalarMatrix TEqn(T, dimless*dimVolume/(dimTime));

    scalar relaxCoeff(0.0);

 Info<< "reading clockProperties\n" << endl;

    IOdictionary clockProperties
    (
            IOobject
            (
                "clockProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
    );

    autoPtr<clockModel> myClock
    (
            clockModel::New
            (
                clockProperties,
                mesh.time()
            )
    );

