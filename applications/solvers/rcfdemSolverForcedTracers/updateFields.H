recurrenceBase.recM().exportVolVectorField("U",URec);

if (useRecP)
{
    recurrenceBase.recM().exportVolScalarField("p",pRec);
}

if (useRecK)
{
    recurrenceBase.recM().exportVolScalarField("k",kRec);
// in case database contains the velocity variance instead of k, do
// kRec *= 0.5;
}
