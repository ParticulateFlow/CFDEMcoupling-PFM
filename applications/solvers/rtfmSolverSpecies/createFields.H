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
    volVectorField U1Rec
    (
        IOobject
        (
            "U.1Rec",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volScalarField alpha1Rec
    (
        IOobject
        (
            "alpha.1Rec",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volVectorField U2Rec
    (
        IOobject
        (
            "U.2Rec",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // calculated fields
    Info<< "\nCreating fields subject to calculation\n" << endl;
    volScalarField alpha2Rec
    (
        IOobject
        (
            "alpha.2Rec",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1-alpha1Rec
    );
        
    
    
    // write fields for t=t_start
    alpha2Rec.write();
//===============================


Info<< "Calculating face flux field phi\n" << endl;
surfaceScalarField phi2Rec
 (
     IOobject
     (
        "phi.2Rec",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
     linearInterpolate(U2Rec*alpha2Rec) & mesh.Sf()
 );
phi2Rec.write();
 
  singlePhaseTransportModel laminarTransport(U2Rec, phi2Rec);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U2Rec, phi2Rec, laminarTransport)
    );
    

    // transport stuff

 // create concentration field
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
    
    T.write();

    scalar relaxCoeff(0.0);