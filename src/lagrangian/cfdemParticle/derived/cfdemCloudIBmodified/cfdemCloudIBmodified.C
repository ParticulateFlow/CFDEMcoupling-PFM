/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-2015 DCS Computing GmbH, Linz
                                Copyright 2015-     JKU Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "fileName.H"
#include "cfdemCloudIBmodified.H"
#include "voidFractionModel.H"
#include "forceModel.H"
#include "locateModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cfdemCloudIBContinuousForcing::cfdemCloudIBContinuousForcing
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    pRefCell_(readLabel(mesh.solutionDict().subDict("PISO").lookup("pRefCell"))),
    pRefValue_(readScalar(mesh.solutionDict().subDict("PISO").lookup("pRefValue"))),
    haveEvolvedOnce_(false),
    skipLagrangeToEulerMapping_(false)
{

    if(this->couplingProperties().found("skipLagrangeToEulerMapping"))
    {
        Info << "Will skip lagrange-to-Euler mapping..." << endl;
        skipLagrangeToEulerMapping_=true;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfdemCloudIBContinuousForcing::~cfdemCloudIBContinuousForcing()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool cfdemCloudIBContinuousForcing::evolve()
{
    numberOfParticlesChanged_ = false;
    arraysReallocated_=false;
    bool doCouple=false;

    if (dataExchangeM().doCoupleNow())
    {
        Info << "\n timeStepFraction() = " << dataExchangeM().timeStepFraction() << endl;
        dataExchangeM().couple(0);
        doCouple=true;

        if(!skipLagrangeToEulerMapping_ || !haveEvolvedOnce_)
        {
          if(verbose_) Info << "- getDEMdata()" << endl;
          getDEMdata();
          Info << "nr particles = " << numberOfParticles() << endl;

          // search cellID of particles
          if(verbose_) Info << "- findCell()" << endl;
          locateM().findCell(NULL,positions_,cellIDs_,numberOfParticles());
          if(verbose_) Info << "findCell done." << endl;

          // set void fraction field
          if(verbose_) Info << "- setvoidFraction()" << endl;
          voidFractionM().setvoidFraction(NULL,voidfractions_,particleWeights_,particleVolumes_,particleV_);
          if(verbose_) Info << "setvoidFraction done." << endl;
        }

        // set particles forces
        if(verbose_) Info << "- setForce(forces_)" << endl;
        for(int index = 0;index <  numberOfParticles_; ++index){
            for(int i=0;i<3;i++){
                impForces_[index][i] = 0;
                expForces_[index][i] = 0;
                DEMForces_[index][i] = 0;
            }
        }
        for (int i=0;i<nrForceModels();i++) forceM(i).setForce();
        if(verbose_) Info << "setForce done." << endl;

        // write DEM data
        if(verbose_) Info << " -giveDEMdata()" << endl;
        giveDEMdata();

        dataExchangeM().couple(1);

        haveEvolvedOnce_=true;
    }
    Info << "evolve done." << endl;

    // do particle IO
    IOM().dumpDEMdata();

    return doCouple;
}

void cfdemCloudIBContinuousForcing::calcForcingTerm(volVectorField& Us)
{
    // set particle velocity field
    if(verbose_) Info << "- setVelocity(velocities_)" << endl;
    label cell = 0;
    vector uP(0,0,0);
    for (int index = 0; index < numberOfParticles(); ++index)
    {
        for (int subCell = 0; subCell < voidFractionM().cellsPerParticle()[index][0]; subCell++)
        {
            cell = cellIDs()[index][subCell];
            if (cell >= 0)
            {
                // calc particle velocity
                for(int i=0;i<3;i++) uP[i] = velocities()[index][i];
                Us[cell] = (1-voidfractions_[index][subCell])*uP;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
