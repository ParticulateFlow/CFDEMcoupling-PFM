if (UMax.value() > 0)
{
    forAll(U,cellI)
    {
        scalar mU(mag(U[cellI]));
        if (mU > UMax.value())
        {
            U[cellI] *= UMax.value() / mU;
        }
    }
}
