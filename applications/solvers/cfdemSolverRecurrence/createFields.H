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
    
    volScalarField voidfractionRec
    (
        IOobject
        (
            "voidfractionRec",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volVectorField UsRec
    (
        IOobject
        (
            "UsRec",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // calculated fields
    Info<< "\nCreating fields subject to calculation\n" << endl;
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

#ifndef createPhi_H
#define createPhi_H
Info<< "Calculating face flux field phi\n" << endl;
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
     linearInterpolate(URec*voidfractionRec) & mesh.Sf()
 );
phiRec.write();
#endif
 
  singlePhaseTransportModel laminarTransport(URec, phiRec);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(URec, phiRec, laminarTransport)
    );
