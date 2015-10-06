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
    
    // calculated fields
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedVector("U",dimVelocity,vector(0,0,0))
    );
    
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

#ifndef createPhi_H
#define createPhi_H
Info<< "Reading/calculating face flux field phi\n" << endl;
// create phi in recModel and put reference here
#endif