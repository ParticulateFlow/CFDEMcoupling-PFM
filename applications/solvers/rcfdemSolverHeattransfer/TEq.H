    volScalarField rhoeps = rhoRec*voidfractionRec;

    particleCloud.energyContributions(Qsource);

    thDiff = particleCloud.thermCondM().thermDiff();

    dimensionedScalar Cp = particleCloud.thermCondM().Cp()*718/1005;

    K = 0.5*magSqr(URec);

    addSource = fvc::ddt(rhoeps, K) + fvc::div(phiRec, K) +
    fvc::div
    (
        fvc::absolute(phiRec/fvc::interpolate(rhoRec), voidfractionRec*URec), pRec
    );

    fvScalarMatrix TEqn =
    (
        fvm::ddt(rhoeps, T)
            + fvm::div(phiRec, T)
            - fvm::laplacian(rhoeps*thDiff, T)
            + addSource/Cp
            - Qsource/Cp
            ==
            fvOptions(rhoeps, T)    // no fvOptions support yet
    );

    //TEqn.relax(relaxCoeff);

    fvOptions.constrain(TEqn); // no fvOptions support yet

    TEqn.solve();
