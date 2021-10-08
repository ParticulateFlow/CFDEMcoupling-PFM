        // get mixture properties
    Cs  = mixture.Cs();
        diffusionCorrection = mixture.diffusionCorrection();
        Deff=particleCloud.diffCoeffM().diffCoeff(); 

        // get scalar source from DEM        
        particleCloud.massContributions(Sm);
        particleCloud.massCoefficients(Smi);

        fvScalarMatrix CEqn
        (
         fvm::ddt(voidfraction,C)
        + fvm::div(phi,C)
         - fvm::laplacian(Deff*voidfraction,C)
        + fvm::div(fvc::interpolate(Deff*voidfraction)*diffusionCorrection*mesh.magSf(), C)
         ==
     Sm + fvm::Sp(Smi,C)
        );

        CEqn.relax();
        fvOptions.constrain(CEqn);
        CEqn.solve();
