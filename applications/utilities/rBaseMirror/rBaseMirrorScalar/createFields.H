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
    word dataBaseName(mirrorProperties.lookup("dataBaseName"));
