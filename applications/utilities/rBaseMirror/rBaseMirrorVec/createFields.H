    IOdictionary mirrorProperties
    (
        IOobject
        (
            "mirrorProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    vector refPoint(mirrorProperties.lookup("refPoint"));
    vector refDirection(mirrorProperties.lookup("refDirection"));

    word fieldName(mirrorProperties.lookup("fieldName"));

    volVectorField origField
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh
    );
