
volScalarField alphaEff("alphaEff", turbulence->nu()/Sc + alphat);

CEqn =
(
    fvm::ddt(C)
    + fvm::div(phiRec, C)
    - fvm::laplacian(alphaEff, C)
    ==
    fvOptions(C)
);

CEqn.relax(relaxCoeff);

fvOptions.constrain(CEqn);

CEqn.solve();

fvOptions.correct(C);
