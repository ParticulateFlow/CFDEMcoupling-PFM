/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
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

#include "cfdemCloudMS.H"
#include "voidFractionModel.H"
#include "forceModelMS.H"
#include "locateModel.H"
#include "dataExchangeModel.H"

//#include "mpi.h" // only for debug reason
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cfdemCloudMS::cfdemCloudMS
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    positionsCM_(NULL),
    velocitiesCM_(NULL),
    cellIDsCM_(NULL),
    bodies_(NULL),
    nrigids_(NULL),
    exCM_(NULL),
    eyCM_(NULL),
    ezCM_(NULL),
    VclumpCM_(NULL),
    SclumpCM_(NULL),
    scalingCM_(NULL),
    typeCM_(NULL),
    Cclump_ex_(NULL),
    Cclump_ey_(NULL),
    impForcesCM_(NULL),
    expForcesCM_(NULL),
    DEMForcesCM_(NULL),
    numberOfClumps_(-1),
    numberOfClumpsChanged_(false),
    forceModels_(couplingProperties_.lookup("forceModelsMS"))
{
    forceModel_ = new autoPtr<forceModelMS>[nrForceModels()];
    for (int i=0;i<nrForceModels();i++)
    {
        forceModel_[i] = forceModelMS::New
        (
            couplingProperties_,
            *this,
            forceModels_[i]
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfdemCloudMS::~cfdemCloudMS()
{
    delete positionsCM_;
    delete velocitiesCM_;
    delete cellIDsCM_;
    delete bodies_;
    delete nrigids_;
    delete exCM_;
    delete eyCM_;
    delete ezCM_;
    delete VclumpCM_;
    delete SclumpCM_;
    delete scalingCM_;
    delete typeCM_;
    delete Cclump_ex_;
    delete Cclump_ey_;
    delete impForcesCM_;
    delete expForcesCM_;
    delete DEMForcesCM_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//             PUBLIC MEMBER FUNCTIONS

void cfdemCloudMS::getDEMdata()
{
    cfdemCloud::getDEMdata();
    // update NClumpTypes in data exch. model
    //dataExchangeM().checkNClumpTypes();

    dataExchangeM().getData("xcm","vector-multisphere",positionsCM_);   // position of the centre of mass
    dataExchangeM().getData("vcm","vector-multisphere",velocitiesCM_);  // velocity of the centre of mass
    dataExchangeM().getData("body","scalar-atom",bodies_);              // clump-particle connex

    dataExchangeM().getData("nrigid","scalar-multisphere",nrigids_);    // # particles in clump
    dataExchangeM().getData("ex_space","vector-multisphere",exCM_);     // axis of inertia
    dataExchangeM().getData("ey_space","vector-multisphere",eyCM_);     // axis of inertia
    dataExchangeM().getData("ez_space","vector-multisphere",ezCM_);     // axis of inertia

//    dataExchangeM().getScalarData("Vclump",VclumpCM_);   // Volume of the clump
//    dataExchangeM().getScalarData("Sclump",SclumpCM_);   // surface area of the clump

//    dataExchangeM().getScalarData("scaling",scalingCM_); // scaling of the clump
//    dataExchangeM().getScalarData("typeCM",typeCM_);    // type of the clump

//    dataExchangeM().getScalarData("Cclump_ex",Cclump_ex_);   // cross section of the clump in ex normal direction
//    dataExchangeM().getScalarData("Cclump_ey",Cclump_ey_);   // cross section of the clump in ey normal direction

    // calc Saequi, surface area of the vol aequivalent sphere
      // Saequi=pi*pow(6*VclumpCM_/pi,2./3.);

    // calc Caequi, cross section of the vol aequivalent sphere
      // Caequi=Saequi/4;
}

void Foam::cfdemCloudMS::giveDEMdata()
{
    for(int index = 0;index < numberOfClumps(); ++index)
    {
        for(int i=0;i<3;i++){
            impForcesCM()[index][i] += expForcesCM()[index][i] + DEMForcesCM()[index][i];
        }
    }
    if(forceM(0).coupleForce()) dataExchangeM().giveData("dragforce","vector-multisphere",impForcesCM());
    if(verbose_) Info << "giveDEMdata done." << endl;
}

bool cfdemCloudMS::evolve
(
    volScalarField& alpha,
    volVectorField& Us,
    volVectorField& U
)
{
    if(cfdemCloud::evolve(alpha,Us,U)) return true;

    return false;
}

bool cfdemCloudMS::reAllocArrays() const
{
    if(cfdemCloud::reAllocArrays())
    {
        int nClumpTypes = dataExchangeM().nClumpTypes();

        // get arrays of new length
        dataExchangeM().allocateArray(positionsCM_,0,3,"nbodies");
        dataExchangeM().allocateArray(velocitiesCM_,0,3,"nbodies");
        dataExchangeM().allocateArray(cellIDsCM_,0,1,"nbodies");
        dataExchangeM().allocateArray(bodies_,0,3);
        dataExchangeM().allocateArray(nrigids_,0,1,"nbodies");
        dataExchangeM().allocateArray(exCM_,0,3,"nbodies");
        dataExchangeM().allocateArray(eyCM_,0,3,"nbodies");
        dataExchangeM().allocateArray(ezCM_,0,3,"nbodies");
        dataExchangeM().allocateArray(VclumpCM_,0,3,nClumpTypes);
        dataExchangeM().allocateArray(SclumpCM_,0,3,nClumpTypes);
        dataExchangeM().allocateArray(scalingCM_,0,3,"nbodies");
        dataExchangeM().allocateArray(typeCM_,0,3,"nbodies");
        dataExchangeM().allocateArray(Cclump_ex_,0,3,nClumpTypes);
        dataExchangeM().allocateArray(Cclump_ey_,0,3,nClumpTypes);
        dataExchangeM().allocateArray(impForcesCM_,0,3,"nbodies");
        dataExchangeM().allocateArray(expForcesCM_,0,3,"nbodies");
        dataExchangeM().allocateArray(DEMForcesCM_,0,3,"nbodies");
        return true;
    }
    return false;
}

void Foam::cfdemCloudMS::setNumberOfParticles(int nP)
{
    cfdemCloud::setNumberOfParticles(nP);

    int nC = dataExchangeM().getNumberOfClumps();

    if(nC != numberOfClumps())
    {
        numberOfClumpsChanged_ = true;
        numberOfClumps_ = nC;
    }
}

void Foam::cfdemCloudMS::findCells()
{
    cfdemCloud::findCells();
    locateM().findCell(NULL,positionsCM_,cellIDsCM_,numberOfClumps());
}

void Foam::cfdemCloudMS::setForces()
{
    cfdemCloud::setForces();

    resetArray(impForces_,numberOfParticles(),3);
    resetArray(expForces_,numberOfParticles(),3);
    resetArray(DEMForces_,numberOfParticles(),3);
    for (int i=0;i<cfdemCloudMS::nrForceModels();i++) cfdemCloudMS::forceM(i).setForce();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//             PUBLIC MEMBER FUNCTIONS

vector Foam::cfdemCloudMS::positionCM(int index)
{
    vector pos;
    for(int i=0;i<3;i++) pos[i] = positionsCM()[index][i];
    return pos;
}

vector Foam::cfdemCloudMS::velocityCM(int index)
{
    vector vel;
    for(int i=0;i<3;i++) vel[i] = velocitiesCM()[index][i];
    return vel;
}

label Foam::cfdemCloudMS::cellIDCM(int index)
{
    return cellIDsCM_[index][0];
}

label Foam::cfdemCloudMS::body(int index)
{
    return bodies_[0][index]-1;
}

label Foam::cfdemCloudMS::nrigid(int index)
{
    return nrigids_[0][index];
//    return nrigids_[0][0];
}

const forceModel& Foam::cfdemCloudMS::forceM(int i)
{
    return forceModel_[i];
}

int Foam::cfdemCloudMS::nrForceModels()
{
    return forceModels_.size();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
