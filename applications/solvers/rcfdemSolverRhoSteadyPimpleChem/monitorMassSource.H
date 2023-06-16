{
    scalar sourceStrength = 0.0;
    forAll(p,cellI)
    {
        sourceStrength += (Sm0[cellI]*p[cellI]+Sm1[cellI])*p.mesh().V()[cellI];
    }
    reduce(sourceStrength, sumOp<scalar>());
    Info << "total mass source strength = " << sourceStrength << endl;
}
