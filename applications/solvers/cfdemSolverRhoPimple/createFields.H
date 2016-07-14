Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<psiThermo> pThermo
    (
        psiThermo::New(mesh)
    );
    psiThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");
    volScalarField& p = thermo.p();
    volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi();

    Info<< "Reading field rho\n" << endl;
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
        thermo.rho()
    );


    
    Info<< "Reading field U\n" << endl;
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

    Info<< "\nCreating fluid-particle heat flux field\n" << endl;
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
    
    Info<< "\nCreating fluid-particle heat flux coefficient field\n" << endl;
    volScalarField Qcoeff
    (
        IOobject
        (
            "Qcoeff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1,-1,-3,-1,0,0,0), 0.0)
    );
    
    Info<< "\nCreating thermal diffusivity field\n" << endl;
    volScalarField thDiff
    (
        IOobject
        (
            "thDiff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,2,-1,0,0,0,0), 0.0)
    );
    
    Info<< "\nCreating body force field\n" << endl;
    volVectorField fOther
    (
        IOobject
        (
            "fOther",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
    );

    #ifndef compressibleCreatePhi_H
    #define compressibleCreatePhi_H
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
         linearInterpolate(rho*U*voidfraction) & mesh.Sf()
    );
    #endif

    dimensionedScalar rhoMax(pimple.dict().lookup("rhoMax"));
    dimensionedScalar rhoMin(pimple.dict().lookup("rhoMin"));

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
    );

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));
    
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
        //dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0), 1.0)
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

//===============================

volScalarField divphi
    (
        IOobject
        (
            "divphi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0), 0.0)
    );
    
    volScalarField divphi2
    (
        IOobject
        (
            "divphi2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0), 0.0)
    );

volVectorField gP
    (
        IOobject
        (
            "gP",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);
    
volVectorField gP2
    (
        IOobject
        (
            "gP2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

    volScalarField rhoeps
    (
        IOobject
        (
            "rhoeps",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(1, -3, 0, 0, 0), 0.0)
    );
