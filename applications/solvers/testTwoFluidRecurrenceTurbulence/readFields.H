// update flow fields
recurrenceBase.recM().exportVolScalarField("alpha."+carrierPhaseName,alpha2);
recurrenceBase.recM().exportVolScalarField("alpha."+dispersePhaseName,alpha1);

recurrenceBase.recM().exportVolVectorField("U."+carrierPhaseName,U2);

// update turbulence models
phase1.correctTurbulence();
phase2.correctTurbulence();
