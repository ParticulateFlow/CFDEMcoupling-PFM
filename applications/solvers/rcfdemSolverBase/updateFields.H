// is it neccessary to extend recurrence path?
if(recurrenceBase.recM().endOfPath())
{
    recurrenceBase.extendPath();
}

for (label i=0;i<numRecVolScalarFieldConnections;i++)
{
    recurrenceBase.recM().exportVolScalarField(recVolScalarFieldConnections[i][0],volScalarFieldsRec[i]);
}

for (label i=0;i<numRecVolVectorFieldConnections;i++)
{
    recurrenceBase.recM().exportVolVectorField(recVolVectorFieldConnections[i][0],volVectorFieldsRec[i]);
}

for (label i=0;i<numRecSurfaceScalarFieldConnections;i++)
{
    recurrenceBase.recM().exportSurfaceScalarField(recSurfaceScalarFieldConnections[i][0],surfaceScalarFieldsRec[i]);
}

