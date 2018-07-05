{
    molarConc = 0.0 * molarConc;
    forAll(Y, i)
    {
        volScalarField& Yi = Y[i];
        dimensionedScalar mi("mi",dimensionSet(1, 0, 0, 0, -1),composition.W(i));
        mi /= 1000.0; // g to kg
        molarConc += rho * Yi / mi;
    }
}

// ************************************************************************* //
