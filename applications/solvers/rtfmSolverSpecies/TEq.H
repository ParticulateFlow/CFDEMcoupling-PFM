    TEqn =
    (
        fvm::ddt(alpha2Rec, T)
        + fvm::div(phi2Rec, T)
        - fvm::laplacian(alpha2Rec*turbulence->nu(), T)
        ==
        fvOptions(alpha2Rec, T) // no fvOptions support yet
    );

    TEqn.relax(relaxCoeff);

    fvOptions.constrain(TEqn); // no fvOptions support yet

    TEqn.solve();
