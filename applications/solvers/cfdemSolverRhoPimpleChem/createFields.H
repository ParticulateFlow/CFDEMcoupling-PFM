    //  thermodynamics, chemistry

#if OPENFOAM_VERSION_MAJOR < 6
    Info<< "Creating combustion model\n" << endl;
    autoPtr<combustionModels::rhoCombustionModel> combustion
    (
        combustionModels::rhoCombustionModel::New(mesh)
    );
    rhoReactionThermo& thermo = combustion->thermo();
#else
    Info<< "Reading thermophysical properties\n" << endl;
    autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
    rhoReactionThermo& thermo = pThermo();
#endif
    thermo.validate(args.executable(), "h", "e");

    basicSpecieMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    // read molecular weight
#if OPENFOAM_VERSION_MAJOR < 6
    volScalarField W(composition.W());
#else
    volScalarField W(thermo.W());
#endif

    Switch propagateInertSpecie(true);

    const word inertSpecie(thermo.lookup("inertSpecie"));

    const scalar inertLowerBound(thermo.lookupOrDefault<scalar>("inertLowerBound",0.0));

    const scalar inertUpperBound(thermo.lookupOrDefault<scalar>("inertUpperBound",1.0));

    if (!composition.contains(inertSpecie))
    {
        FatalErrorIn(args.executable())
            << "Specified inert specie '" << inertSpecie << "' not found in "
            << "species list. Available species:" << composition.species()
            << exit(FatalError);
    }

    Info<< "inert will be bounded in [" << inertLowerBound << "," << inertUpperBound << "]" << endl;

    volScalarField& p = thermo.p();

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(thermo.he());

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    //  kinematic fields
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

    volScalarField rhoeps ("rhoeps", rho*voidfraction);


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
    volScalarField QCoeff
    (
        IOobject
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

    Info<< "\nCreating thermal conductivity field\n" << endl;
    volScalarField thCond
    (
        IOobject
        (
            "thCond",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1,1,-3,-1,0,0,0), 0.0),
        "zeroGradient"
    );

    Info<< "\nCreating heat capacity field\n" << endl;
    volScalarField Cpv
    (
        IOobject
        (
            "Cpv",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,2,-2,-1,0,0,0), 0.0)
    );

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

    dimensionedScalar rhoMax
    (
        dimensionedScalar::lookupOrDefault
        (
            "rhoMax",
            pimple.dict(),
            dimDensity,
            GREAT
        )
    );

    dimensionedScalar rhoMin
    (
        dimensionedScalar::lookupOrDefault
        (
            "rhoMin",
            pimple.dict(),
            dimDensity,
            0
        )
    );

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

#if OPENFOAM_VERSION_MAJOR >= 6
    Info<< "Creating combustion model\n" << endl;
    autoPtr<CombustionModel<rhoReactionThermo>> combustion
    (
        CombustionModel<rhoReactionThermo>::New(thermo, turbulence())
    );
#endif

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

#if OPENFOAM_VERSION_MAJOR < 5
    volScalarField dQ
    (
        IOobject
        (
            "dQ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
    );
#else
    volScalarField Qdot
    (
        IOobject
        (
            "Qdot",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
    );
#endif

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

    volScalarField molarConc
    (
        IOobject
        (
            "molarConc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0, -3, 0, 0, 1),0)
    );

    volScalarField dSauter
    (
        IOobject
        (
            "dSauter",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0, 1, 0, 0, 0,0,0),0)
    );
//===============================
