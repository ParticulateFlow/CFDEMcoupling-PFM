{
    scalar massSink = 0.0;
    forAll(Yi,cellI)
    {
        if (Yi[cellI] <= 0.0)
        {
            massSink += rhoeps[cellI]*Yi[cellI]*Yi.mesh().V()[cellI];
        }
    }
    reduce(massSink, sumOp<scalar>());
    Info << Y[i].name() << ": mass sink = " << massSink << endl;
}
