// get current weights for various databases

// A: triggered over current value of boundary field
// word boundaryName = "inlet";
// label myinlet = mesh.boundary().findPatchID(boundaryName);
// label startIndex = mesh.boundary()[boundaryName].start();

// B: explicitly define weights

scalarList wList(weightDict.lookupOrDefault("weights",scalarList(1,0.0)));

recurrenceBase.recP().updateIntervalWeights(wList);

// is it neccessary to extend recurrence path?
if(recurrenceBase.recM().endOfPath())
{
    recurrenceBase.extendPath();
}

recurrenceBase.recM().exportVolScalarField("voidfraction",voidfractionRec);
recurrenceBase.recM().exportVolScalarField("p",pRec);
recurrenceBase.recM().exportVolVectorField("Us",UsRec);
recurrenceBase.recM().exportSurfaceScalarField("phi",phiRec);

Info << "current database weights: = " << wList << endl;
Info << "current database: " << recurrenceBase.recM().currDataBase() << endl;
for(int i=0;i<wList.size();i++)
{
    scalar w = wList[i];
    if (recurrenceBase.recM().currDataBase() == i) w -= 1.0;
    phiRec += w*recurrenceBase.recM().exportSurfaceScalarFieldAve("phi",i)();
}

{
    volScalarField& NuField(const_cast<volScalarField&>(mesh.lookupObject<volScalarField> ("NuField")));
    recurrenceBase.recM().exportVolScalarField("NuField",NuField);
}

