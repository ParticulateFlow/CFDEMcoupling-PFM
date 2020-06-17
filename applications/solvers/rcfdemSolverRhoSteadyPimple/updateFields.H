// is it neccessary to extend recurrence path?
if(recurrenceBase.recM().endOfPath())
{
    recurrenceBase.extendPath();
}

recurrenceBase.recM().exportVolScalarField("voidfraction",voidfractionRec);
recurrenceBase.recM().exportVolVectorField("Us",UsRec);
