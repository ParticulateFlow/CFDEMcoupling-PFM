    IOdictionary recProperties
    (
        IOobject
        (
            "recProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    scalar threshold(readScalar(recProperties.lookup("threshold")));

    volVectorField U_smooth
    (
        IOobject
        (
            "U_smooth",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0,0,0), vector(0.0, 0.0, 0.0))
    );
