    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading physical velocity field U" << endl;
    Info<< "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//========================
// drag law modelling
//========================

    Info<< "\nReading momentum exchange field Ksl\n" << endl;
    volScalarField Ksl
    (
        IOobject
        (
            "Ksl",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "\nReading voidfraction field voidfraction = (Vgas/Vparticle)\n" << endl;
    volScalarField voidfraction
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "\nCreating density field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading particle velocity field Us\n" << endl;
    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//========================
// scalar field modelling
//========================

    Info<< "\nCreating temperature field\n" << endl;
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

    Info<< "\nCreating fluid-particle heat flux field\n" << endl;
    volScalarField Qsource
    (
        IOobject
        (
            "Qsource",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    Info<< "\nCreating thermal diffusivity field\n" << endl;
    volScalarField thDiff
    (
        IOobject
        (
            "thDiff",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    Info<< "\nCreating temperature gradient field\n" << endl;
    volVectorField gradT
    (
        IOobject
        (
            "gradT",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
    dimensionedScalar Cp(transportProperties.lookup("Cp"));


//========================

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
     linearInterpolate(U*voidfraction) & mesh.Sf()
 );
#endif



    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );
