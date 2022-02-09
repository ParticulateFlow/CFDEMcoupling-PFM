    // dummy fields
    Info << "\nCreating dummy pressure and density fields\n" << endl;
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
    Info << "\nCreating recurrence fields.\n" << endl;
    volVectorField URec
    (
        IOobject
        (
            "URec",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("URec", dimensionSet(0, 1, -1, 0, 0), vector::zero)
    );

    Switch updateURec(false);
    if (URec.headerOk())
    {
        updateURec = true;
        URec.writeOpt() = IOobject::AUTO_WRITE;
    }

    volScalarField voidfractionRec
    (
        IOobject
        (
            "voidfractionRec",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("voidfractionRec", dimensionSet(0, 0, 0, 0, 0), 1.0)
    );

    Switch updateVoidfractionRec(false);
    if (voidfractionRec.headerOk())
    {
        updateVoidfractionRec = true;
        voidfractionRec.writeOpt() = IOobject::AUTO_WRITE;
    }

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
        dimensionedVector("URec", dimensionSet(0, 1, -1, 0, 0), vector::zero)
    );

    Switch updateUsRec(false);
    if (UsRec.headerOk())
    {
        updateUsRec = true;
        UsRec.writeOpt() = IOobject::AUTO_WRITE;
    }

    // calculated fields
    Info << "\nCreating fields subject to calculation\n" << endl;
    volScalarField voidfraction
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        UsRec
    );

    // write fields for t=t_start
    voidfraction.write();
    Us.write();
//===============================

    Info << "Calculating face flux field phi\n" << endl;
    surfaceScalarField phiRec
    (
        IOobject
        (
            "phiRec",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        linearInterpolate(URec*voidfractionRec) & mesh.Sf()
    );

    Switch updatePhiRec(false);
    if (phiRec.headerOk())
    {
        updatePhiRec = true;
        phiRec.writeOpt() = IOobject::AUTO_WRITE;
        phiRec.write();
    }

    singlePhaseTransportModel laminarTransport(URec, phiRec);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(URec, phiRec, laminarTransport)
    );

    IOdictionary recDict
    (
        IOobject
        (
            "recProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word voidfractionFieldName(recDict.lookupOrDefault<word>("voidfractionFieldName","voidfraction"));
    word UFieldName(recDict.lookupOrDefault<word>("UFieldName","U"));
    word UsFieldName(recDict.lookupOrDefault<word>("UsFieldName","Us"));
    word fluxFieldName(recDict.lookupOrDefault<word>("fluxFieldName","phi"));


   // place to put weight functions
    IOdictionary weightDict
    (
        IOobject
        (
            "weightDict",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    if (!weightDict.headerOk())
    {
        weightDict.add("weights",scalarList(1,1.0));
    }
    scalarList weights(weightDict.lookup("weights"));
    Info << "database initial weights: " << weights << endl;
