scalarList wList(weightDict.lookupOrDefault("weights",scalarList(1,1.0)));

recurrenceBase.recP().updateIntervalWeights(wList);

if(recurrenceBase.recM().endOfPath())
{
    recurrenceBase.extendPath();
}

// update fields where necessary
if (updateVoidfractionRec)
{
    recurrenceBase.recM().exportVolScalarField(voidfractionFieldName,voidfractionRec);
}

if (updateURec)
{
    recurrenceBase.recM().exportVolVectorField(UFieldName,URec);
}

if (updateUsRec)
{
    recurrenceBase.recM().exportVolVectorField(UsFieldName,UsRec);
}

if (updatePhiRec)
{
    recurrenceBase.recM().exportSurfaceScalarField(fluxFieldName,phiRec);
}
