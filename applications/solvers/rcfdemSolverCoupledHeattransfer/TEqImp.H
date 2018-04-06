    volScalarField rhoeps = rhoRec*voidfractionRec;

    particleCloud.energyContributions(Qsource);

    particleCloud.energyCoefficients(QCoeff);

    thCond = particleCloud.thermCondM().thermCond();

    dimensionedScalar Cv = particleCloud.thermCondM().Cv();

    //K = 0.5*magSqr(URec);

    addSource = fvc::div(phiRec/fvc::interpolate(rhoRec), pRec);
    // main contribution due to gas expansion, not due to transport of kinetic energy
    // fvc::ddt(rhoeps, K) + fvc::div(phiRec, K)

    fvScalarMatrix TEqn =
    (
        fvm::ddt(rhoeps, T)
        + fvm::div(phiRec, T)
        + addSource/Cv
        - fvm::laplacian(voidfractionRec*thCond/Cv, T)
        - Qsource/Cv
        - fvm::Sp(QCoeff/Cv, T)
        ==
        fvOptions(rhoeps, T)    // no fvOptions support yet
    );

    fvOptions.constrain(TEqn); // no fvOptions support yet

    TEqn.solve();

    particleCloud.clockM().start(31,"postFlow");
    counter++;

    if((counter - couplingSubStep) % dtDEM2dtCFD == 0)
        particleCloud.postFlow();

    particleCloud.clockM().stop("postFlow");
