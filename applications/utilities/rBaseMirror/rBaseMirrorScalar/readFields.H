
    volScalarField origField
    (
        IOobject
        (
            fieldName,
            recTime.timePath(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
