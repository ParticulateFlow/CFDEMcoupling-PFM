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
