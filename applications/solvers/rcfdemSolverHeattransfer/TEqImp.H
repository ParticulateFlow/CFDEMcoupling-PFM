    volScalarField rhoeps = rhoRec*voidfractionRec;

    particleCloud.energyContributions(Qsource);

    particleCloud.energyCoefficients(QCoeff);

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
        + addSource/Cv
        - fvm::laplacian(voidfractionRec*thCond/Cv, T)
        - Qsource/Cv
        - fvm::Sp(QCoeff/Cv, T)
        ==
        fvOptions(rhoeps, T)    // no fvOptions support yet
    );

    //TEqn.relax(relaxCoeff);

    fvOptions.constrain(TEqn); // no fvOptions support yet

    TEqn.solve();

    particleCloud.clockM().start(31,"postFlow");
    counter++;

    if((counter - couplingSubStep) % dtDEM2dtCFD == 0)
        particleCloud.postFlow();

    particleCloud.clockM().stop("postFlow");
