// case 1

forAll(rhoRec,cellI)
{
    if (mesh.C()[cellI].z() < 0.00228)
        rhoRec[cellI] = 1.18+(1.085-1.18)*Foam::exp(-0.065*runTime.timeOutputValue());
    else if (mesh.C()[cellI].z() < 0.00456)
        rhoRec[cellI] = 1.18+(1.01-1.18)*Foam::exp(-0.05*runTime.timeOutputValue());
    else if (mesh.C()[cellI].z() < 0.00684)
        rhoRec[cellI] = 1.18+(0.98-1.18)*Foam::exp(-0.0425*runTime.timeOutputValue());
    else
        rhoRec[cellI] = 1.18+(0.955-1.18)*Foam::exp(-0.0425*runTime.timeOutputValue());   
}


// case 2
/*
forAll(rhoRec,cellI)
{
    if (mesh.C()[cellI].z() < 0.00228)
        rhoRec[cellI] = 1.18+(1.115-1.18)*Foam::exp(-0.065*runTime.timeOutputValue());
    else if (mesh.C()[cellI].z() < 0.00456)
        rhoRec[cellI] = 1.18+(1.04-1.18)*Foam::exp(-0.05*runTime.timeOutputValue());
    else if (mesh.C()[cellI].z() < 0.00684)
        rhoRec[cellI] = 1.18+(1.005-1.18)*Foam::exp(-0.0425*runTime.timeOutputValue());
    else
        rhoRec[cellI] = 1.18+(0.96-1.18)*Foam::exp(-0.0425*runTime.timeOutputValue());   
}
*/