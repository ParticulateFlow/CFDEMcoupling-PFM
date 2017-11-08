recurrenceBase.recM().exportVolScalarField("alpha.air",alpha1Rec);
alpha2Rec=1-alpha1Rec;
recurrenceBase.recM().exportVolVectorField("U.air",U1Rec);
recurrenceBase.recM().exportVolVectorField("U.water",U2Rec);
recurrenceBase.recM().exportSurfaceScalarField("phi.water",phi2Rec);

