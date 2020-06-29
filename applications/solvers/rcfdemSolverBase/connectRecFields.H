// look for connection between fields in database and fields used in solver
// corresponding entry in recProperties:
// recVolScalarFieldConnections ((voidfraction voidfractionRec) (p pRec) ... )

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

List<wordList> recVolScalarFieldConnections(recProperties.lookup("recVolScalarFieldConnections"));
List<wordList> recVolVectorFieldConnections(recProperties.lookup("recVolVectorFieldConnections"));
List<wordList> recSurfaceScalarFieldConnections(recProperties.lookup("recSurfaceScalarFieldConnections"));

label numRecVolScalarFieldConnections = recVolScalarFieldConnections.size();
label numRecVolVectorFieldConnections = recVolVectorFieldConnections.size();
label numRecSurfaceScalarFieldConnections = recSurfaceScalarFieldConnections.size();

PtrList<volScalarField> volScalarFieldsRec;
PtrList<volVectorField> volVectorFieldsRec;
PtrList<surfaceScalarField> surfaceScalarFieldsRec;

for (label i=0;i<numRecVolScalarFieldConnections;i++)
{
    volScalarField &refField(mesh.lookupObjectRef<volScalarField>(recVolScalarFieldConnections[i][1]));
    volScalarFieldsRec.append(&refField);
    Info << "connected " << recVolScalarFieldConnections[i][0] << " and " << recVolScalarFieldConnections[i][1] << endl;

}

for (label i=0;i<numRecVolVectorFieldConnections;i++)
{
    volVectorField &refField(mesh.lookupObjectRef<volVectorField>(recVolVectorFieldConnections[i][1]));
    volVectorFieldsRec.append(&refField);
    Info << "connected " << recVolVectorFieldConnections[i][0] << " and " << recVolVectorFieldConnections[i][1] << endl;
}

for (label i=0;i<numRecSurfaceScalarFieldConnections;i++)
{
    surfaceScalarField &refField(mesh.lookupObjectRef<surfaceScalarField>(recSurfaceScalarFieldConnections[i][1]));
    surfaceScalarFieldsRec.append(&refField);
    Info << "connected " << recSurfaceScalarFieldConnections[i][0] << " and " << recSurfaceScalarFieldConnections[i][1] << endl;

}
