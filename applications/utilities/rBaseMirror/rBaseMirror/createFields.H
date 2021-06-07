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

    word dataBaseName(mirrorProperties.lookup("dataBaseName"));

    const wordList volScalarFieldNames(mirrorProperties.lookup("volScalarFields"));
    const wordList volVectorFieldNames(mirrorProperties.lookup("volVectorFields"));
