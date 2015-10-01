    // dummy fields
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
        dimensionedScalar("rho", dimensionSet(1, -3, 0, 0, 0), 1.0)
    );

    // recurrence fields
    volVectorField URec
    (
        IOobject
        (
            "URec",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedVector("URec",dimVelocity,vector(0,0,0))
    );

    volScalarField voidfractionRec
    (
        IOobject
        (
            "voidfractionRec",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar("voidfractionRec", dimensionSet(0, -3, 0, 0, 0), 1.0)
    );

    volVectorField UsRec
    (
        IOobject
        (
            "UsRec",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedVector("UsRec",dimVelocity,vector(0,0,0))
    );
    
    // calculated fields
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
        dimensionedScalar("voidfraction", dimensionSet(0, -3, 0, 0, 0), 1.0)
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
        dimensionedVector("Us",dimVelocity,vector(0,0,0))
    );

//===============================

//#   include "createPhi.H"
#ifndef createPhi_H
#define createPhi_H
Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
 (
     IOobject
     (
        "phi",
        runTime.timeName(),
         mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
     linearInterpolate(URec*voidfractionRec) & mesh.Sf()
 );
#endif