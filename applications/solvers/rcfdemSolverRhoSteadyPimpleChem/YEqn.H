particleCloud.clockM().start(29,"Y");

tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    combustion->correct();
#if OPENFOAM_VERSION_MAJOR < 5
    dQ = combustion->dQ();
#else
    Qdot = combustion->Qdot();
#endif
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (Y[i].name() == inertSpecie) inertIndex = i;
        if (Y[i].name() != inertSpecie || propagateInertSpecie)
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(voidfractionRec*turbulence->muEff(), Yi)
              ==
                combustion->R(Yi)
              + particleCloud.chemistryM(0).Smi(i)*p/p.prevIter()
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            Yi.relax();

            fvOptions.correct(Yi);

            Yi.max(0.0);
            if (Y[i].name() != inertSpecie) Yt += Yi;
        }
    }

    if (inertIndex!=-1)
    {
        Y[inertIndex].max(inertLowerBound);
        Y[inertIndex].min(inertUpperBound);
    }

    if (propagateInertSpecie)
    {
        if (inertIndex!=-1) Yt /= (1-Y[inertIndex] + ROOTVSMALL);
        forAll(Y,i)
        {
            if (i!=inertIndex)
            {
                volScalarField& Yi = Y[i];
                Yi = Yi/(Yt+ROOTVSMALL);
            }
        }
    }
    else
    {
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].max(0.0);
    }
}

particleCloud.clockM().stop("Y");
