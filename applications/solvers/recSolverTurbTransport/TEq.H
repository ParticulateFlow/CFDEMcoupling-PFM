
volScalarField alphaEff("alphaEff", turbulence->nu()/Sc + dU2/Sct);

TEqn =
(
    fvm::ddt(T)
    + fvm::div(phiRec, T)
    - fvm::laplacian(alphaEff, T)
    ==
    fvOptions(T)
);

TEqn.relax(relaxCoeff);

fvOptions.constrain(TEqn);

TEqn.solve();
