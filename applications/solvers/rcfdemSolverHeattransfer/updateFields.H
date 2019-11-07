recurrenceBase.recM().exportVolScalarField("voidfraction",voidfractionRec);
recurrenceBase.recM().exportVolScalarField("rho",rhoRec);
recurrenceBase.recM().exportVolScalarField("p",pRec);
recurrenceBase.recM().exportVolVectorField("Us",UsRec);
recurrenceBase.recM().exportVolVectorField("U",URec);
recurrenceBase.recM().exportSurfaceScalarField("phi",phiRec);

{
    volScalarField& NuField(const_cast<volScalarField&>(mesh.lookupObject<volScalarField> ("NuField")));
    recurrenceBase.recM().exportVolScalarField("NuField",NuField);
}


#include "updateRho.H"