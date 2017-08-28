{
            volScalarField artMass = rhoeps;
            scalar lowestValue(0.0);
            label lVCell(-1);
            forAll(Yi,cellI)
            {
                if(Yi[cellI] < 0.0)
                {
                   artMass[cellI] *= Yi[cellI];
                   if(artMass[cellI] < lowestValue)
                   {
                       lowestValue=artMass[cellI];
                       lVCell = cellI;
                   }
                }
                else
                {
                    artMass[cellI] *=0.0;
                }
            }
            Info << "\nartificial mass of " << Y[i].name() << " per time step: "<< fvc::domainIntegrate(artMass) << endl;
            if(lVCell > -1)
            {
                Info << "lowest value of " << lowestValue << " at cell " << lVCell << "with coordinates" << endl;
                Info << "\t" << mesh.C()[lVCell].component(0) << "\t" << mesh.C()[lVCell].component(1)  << "\t" << mesh.C()[lVCell].component(2) << endl;
            }
}
