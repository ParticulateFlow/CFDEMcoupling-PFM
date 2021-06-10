scalarList wList(weightDict.lookupOrDefault("weights",scalarList(1,1.0)));

recurrenceBase.recP().updateIntervalWeights(wList);

// is it neccessary to extend recurrence path?
if(recurrenceBase.recM().endOfPath())
{
    recurrenceBase.extendPath();
}

// update fields where necessary
if (updateVoidfractionRec)
{
    recurrenceBase.recM().exportVolScalarField("voidfraction",voidfractionRec);
}

if (updateURec)
{
    recurrenceBase.recM().exportVolVectorField("U",URec);
}

if (updateUsRec)
{
    recurrenceBase.recM().exportVolVectorField("Us",UsRec);
}

if (updatePhiRec)
{
    recurrenceBase.recM().exportSurfaceScalarField("phi",phiRec);
}
