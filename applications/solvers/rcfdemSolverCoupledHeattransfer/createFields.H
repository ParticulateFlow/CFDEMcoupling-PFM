    // dummy fields
    Info << "\nCreating dummy pressure field\n" << endl;
    volScalarField pRec
    (
        IOobject
        (
            "pRec",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    );

    // recurrence fields
    Info << "\nCreating recurrence fields.\n" << endl;

    volScalarField rhoRec
    (
        IOobject
        (
            "rhoRec",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0), 1.0)
    );

    volVectorField URec
    (
        IOobject
        (
            "URec",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0), vector::zero)
    );

    volScalarField voidfractionRec
    (
        IOobject
        (
            "voidfractionRec",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    );

    volVectorField UsRec
    (
        IOobject
        (
            "UsRec",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0), vector::zero)
    );

    // heat transfer fields
    Info << "\nCreating heat transfer fields.\n" << endl;

    volScalarField Qsource
    (
        IOobject
        (
            "Qsource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    );

    volScalarField QCoeff
    (   IOobject
        (
            "QCoeff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1,-1,-3,-1,0,0,0), 0.0)
    );

    volScalarField thCond
    (
        IOobject
        (
            "thCond",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1,1,-3,-1,0,0,0), 0.0)
    );

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

    // calculated fields
    Info << "\nCreating fields subject to calculation\n" << endl;
    volScalarField voidfraction
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        voidfractionRec
    );

    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        UsRec
    );

    // write fields for t=t_start
    voidfraction.write();
    Us.write();
//===============================


    Info << "Calculating face flux field phiRec\n" << endl;
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
        mesh,
        dimensionedScalar("zero", dimensionSet(1,0,-1,0,0,0,0), 0.0)
    );
    phiRec.write();


    Info << "Creating dummy turbulence model\n" << endl;
    singlePhaseTransportModel laminarTransport(URec, phiRec);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(URec, phiRec, laminarTransport)
    );


    const IOdictionary& transportProps = mesh.lookupObject<IOdictionary>("transportProperties");
    dimensionedScalar molMass(transportProps.lookup("molM"));
    // need to scale R down with 1e3 because return value of RR in g, not kg
    dimensionedScalar R("R",dimensionSet(0,2,-2,-1,0,0,0),Foam::constant::thermodynamic::RR / (1e3*molMass.value()));
    Info << "specific gas constant R = " << R << endl;


    volScalarField addSource
    (
        IOobject
        (
           "addSource",
           runTime.timeName(),
           mesh,
           IOobject::READ_IF_PRESENT,
           IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    );


   // place to put weight functions
    IOdictionary weightDict
    (
        IOobject
        (
            "weightDict",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );
    weightDict.add("weights",scalarList(1,1.0));
