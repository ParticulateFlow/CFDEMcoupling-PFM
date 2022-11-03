
currTimeIndex = recBases[0].recM().currentTimeIndex();

Info << "current Time Index = " << currTimeIndex << endl;

recBases[0].recM().exportVolVectorField(Ucoh_pair,UcohRec);
recBases[0].recM().exportVolScalarField(kSGScoh_pair,kcohRec);


label incTimeIndex = incPairTimeIndex_[currTimeIndex][1];

Info << "  incoherent pair Time Index = " << incTimeIndex  << endl;

UincRec = recBases[1].recM().exportVolVectorField(Uinc_pair,incTimeIndex);
kincRec = recBases[1].recM().exportVolScalarField(kSGSinc_pair,incTimeIndex);

kRec = kcohRec+kincRec;

forAll(kRec, cellI)
{
    if (kRec[cellI] < SMALL)
    {
        kRec[cellI] = 0.0;
    }
}

const fvPatchList& patches = mesh.boundary();
forAll(patches, patchI)
{
        kRec.boundaryFieldRef()[patchI] = 0.0;
}

URec = UcohRec + UincRec;
phiRec = linearInterpolate(URec) & mesh.Sf();

nutRec = sqrt(kRec)*delta*0.094;

alphat = nutRec/Sct;
alphat.correctBoundaryConditions();

