{
    scalar sourceStrength = 0.0;
    forAll(Yi,cellI)
    {
        sourceStrength += (Smi0[cellI]*Yi[cellI]+Smi1[cellI])*Yi.mesh().V()[cellI];
    }
    reduce(sourceStrength, sumOp<scalar>());
    Info << Y[i].name() << ": source strength = " << sourceStrength << endl;
}
