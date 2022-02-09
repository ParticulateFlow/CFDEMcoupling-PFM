 // get scalar source from DEM
        particleCloud.forceM(1).manipulateScalarField(Tsource);
        Tsource.correctBoundaryConditions();

        // solve scalar transport equation
        fvScalarMatrix TEqn
        (
           fvm::ddt(voidfraction,T) - fvm::Sp(fvc::ddt(voidfraction),T)
         + fvm::div(phi, T) - fvm::Sp(fvc::div(phi),T)
         - fvm::laplacian(DT*voidfraction, T)
         ==
           Tsource
        );
        TEqn.relax();
        TEqn.solve();
