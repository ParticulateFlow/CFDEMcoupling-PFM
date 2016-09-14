recurrenceBase.recM().exportVolScalarField("voidfraction",alpha1Rec);
alpha2Rec=1-alpha1Rec;
recurrenceBase.recM().exportVolVectorField("U",U1Rec);
recurrenceBase.recM().exportVolVectorField("Us",U2Rec);
recurrenceBase.recM().exportSurfaceScalarField("phiS",phi2Rec);

