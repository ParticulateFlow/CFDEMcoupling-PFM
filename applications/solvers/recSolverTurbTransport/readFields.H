
recurrenceBase.recM().exportVolScalarField("U2Mean",U2Rec);
recurrenceBase.recM().exportVolVectorField("UMean",URec);
phiRec=linearInterpolate(URec) & mesh.Sf();

dU2=sqrt(0.5*mag(U2Rec - magSqr(URec)))*delta*0.094;

forAll(dU2, cellI)
{
    if (U2Rec[cellI]-magSqr(URec[cellI]) < 0.0)
    {
        dU2[cellI] = 0.0;
    }
}
