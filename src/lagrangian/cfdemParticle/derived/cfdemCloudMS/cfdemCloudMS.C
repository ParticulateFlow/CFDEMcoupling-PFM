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

//#include <mpi.h> // only for debug reason
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
    clumpType_(NULL),
    nClumpTypes_(1),
    clumpVol_(NULL),
    clumpDH_(NULL),
    clumpWeights_(NULL),
    //exCM_(NULL),
    //eyCM_(NULL),
    //ezCM_(NULL),
    //SclumpCM_(NULL),
    //scalingCM_(NULL),
    //Cclump_ex_(NULL),
    //Cclump_ey_(NULL),
    impForcesCM_(NULL),
    expForcesCM_(NULL),
    DEMForcesCM_(NULL),
    numberOfClumps_(-1),
    numberOfClumpsChanged_(false),
    manDHdev_(false),
    dHbyV_(scalarList(0)),
    useforcePerClump_(false),
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
    dataExchangeM().destroy(positionsCM_,3);
    dataExchangeM().destroy(velocitiesCM_,3);
    dataExchangeM().destroy(cellIDsCM_,1);
    dataExchangeM().destroy(bodies_,1);
    dataExchangeM().destroy(nrigids_,1);
    dataExchangeM().destroy(clumpType_,1);
    dataExchangeM().destroy(clumpVol_,1);
    dataExchangeM().destroy(clumpDH_,1);
    dataExchangeM().destroy(clumpWeights_,1);
    dataExchangeM().destroy(impForcesCM_,3);
    dataExchangeM().destroy(expForcesCM_,3);
    dataExchangeM().destroy(DEMForcesCM_,3);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//             PUBLIC MEMBER FUNCTIONS

void cfdemCloudMS::getDEMdata()
{
    cfdemCloud::getDEMdata();

    dataExchangeM().getData("xcm","vector-multisphere",positionsCM_);   // position of the centre of mass
    dataExchangeM().getData("vcm","vector-multisphere",velocitiesCM_);  // velocity of the centre of mass
    dataExchangeM().getData("body","scalar-atom",bodies_);              // clump-particle connex
    dataExchangeM().getData("nrigid","scalar-multisphere",nrigids_);    // # particles in clump
    dataExchangeM().getData("clumptype","scalar-multisphere",clumpType_);   // type of the clump

    nClumpTypes_=dataExchangeM().getNumberOfTypes();                        // nr of clump types
    double* typeVol_=dataExchangeM().getTypeVol();                          // volume of the clump

    //- save clump volume and mass
    double **typeDH(NULL);
    dataExchangeM().allocateArray(typeDH,-1.,1,nClumpTypes()+1);
    if(manDHdev_) // use manually defined dH
    {
        for(int k = 1;k <= nClumpTypes(); k++)
            typeDH[k][0]=dHbyV_[k-1]*typeVol_[k];
    }
    else // calc dH from volAeqivalent shpere
    {
        for(int k = 1;k <= nClumpTypes(); k++)
            typeDH[k][0]=pow(typeVol_[k]*1.9099,1./3.); // 6/pi=1.9099 // calc a hydraulic diameter as d of vol equal sphere
    }

    int ct(0);
    for(int ind = 0;ind < numberOfClumps(); ind++)
    {
        ct=clumpType()[0][ind];
        clumpVol_[ind][0] = typeVol_[ct];
        clumpDH_[ind][0]=typeDH[ct][0];
        //Info << "ct=" << ct << endl;
        //Info << "clumpVol()[ind][0]=" << clumpVol()[ind][0] << endl;
        //Info << "clumpDH()[ind][0]=" << clumpDH()[ind][0] << endl;
    }
    dataExchangeM().destroy(typeDH,1);
    // --

    //dataExchangeM().getData("ex_space","vector-multisphere",exCM_);     // axis of inertia
    //dataExchangeM().getData("ey_space","vector-multisphere",eyCM_);     // axis of inertia
    //dataExchangeM().getData("ez_space","vector-multisphere",ezCM_);     // axis of inertia
    //dataExchangeM().getScalarData("Sclump",SclumpCM_);   // surface area of the clump
    //dataExchangeM().getScalarData("scaling",scalingCM_); // scaling of the clump
    //dataExchangeM().getScalarData("Cclump_ex",Cclump_ex_);   // cross section of the clump in ex normal direction
    //dataExchangeM().getScalarData("Cclump_ey",Cclump_ey_);   // cross section of the clump in ey normal direction
}

void cfdemCloudMS::giveDEMdata()
{
    /*for(int index = 0;index < numberOfClumps(); ++index)
    {
        for(int i=0;i<3;i++){
            impForcesCM()[index][i] += expForcesCM()[index][i] + DEMForcesCM()[index][i];
        }
    }
    if(forceM(0).coupleForce()) dataExchangeM().giveData("dragforce","vector-multisphere",impForcesCM());*/
    if(forceM(0).coupleForce()) dataExchangeM().giveData("dragforce_cm","vector-multisphere",DEMForcesCM());
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

bool cfdemCloudMS::reAllocArrays()
{
    if(cfdemCloud::reAllocArrays())
    {
        // get arrays of new length
        dataExchangeM().allocateArray(positionsCM_,0.,3,"nbodies");
        dataExchangeM().allocateArray(velocitiesCM_,0.,3,"nbodies");
        dataExchangeM().allocateArray(cellIDsCM_,-1,1,"nbodies");
        dataExchangeM().allocateArray(bodies_,0,1);
        dataExchangeM().allocateArray(nrigids_,0,1,"nbodies");
        dataExchangeM().allocateArray(clumpType_,0,1,"nbodies");
        dataExchangeM().allocateArray(clumpVol_,0.,1,"nbodies");
        dataExchangeM().allocateArray(clumpDH_,1.,1,"nbodies");
        dataExchangeM().allocateArray(clumpWeights_,1.,1,"nbodies");
        dataExchangeM().allocateArray(impForcesCM_,0.,3,"nbodies");
        dataExchangeM().allocateArray(expForcesCM_,0.,3,"nbodies");
        dataExchangeM().allocateArray(DEMForcesCM_,0.,3,"nbodies");
        return true;
    }
    return false;
}

void cfdemCloudMS::setNumberOfParticles(int nP)
{
    cfdemCloud::setNumberOfParticles(nP);
    int nC = dataExchangeM().getNumberOfClumps();

    if(nC != numberOfClumps())
    {
        numberOfClumpsChanged_ = true;
        numberOfClumps_ = nC;
    }

    // in case last particle has left an ma-_tag_ms is not up to date
    numberOfClumps_ = min(numberOfParticles(),numberOfClumps_);
}

void cfdemCloudMS::findCells()
{
    cfdemCloud::findCells();
    locateM().findCell(NULL,positionsCM_,cellIDsCM_,numberOfClumps());
}

void cfdemCloudMS::setForces()
{
    resetArray(impForces_,numberOfParticles(),3);
    resetArray(expForces_,numberOfParticles(),3);
    resetArray(DEMForces_,numberOfParticles(),3);

    cfdemCloud::setForces();

    resetArray(impForcesCM_,numberOfClumps(),3);
    resetArray(expForcesCM_,numberOfClumps(),3);
    resetArray(DEMForcesCM_,numberOfClumps(),3);
    for (int i=0;i<cfdemCloudMS::nrForceModels();i++) cfdemCloudMS::forceM(i).setForce();
}

void cfdemCloudMS::setParticleForceField()
{
    // set forces per particle
    cfdemCloud::setParticleForceField();

    // set forces per clump
    // seems to be very unstable - exchange field is very inhomogeneous
    if(useforcePerClump_)
    {
        averagingM().setVectorSumSimple
        (
        forceM(0).impParticleForces(),
        impForcesCM_,
        clumpWeights_,
        numberOfClumps()
        );
        averagingM().setVectorSumSimple
        (
        forceM(0).expParticleForces(),
        expForcesCM_,
        clumpWeights_,
        numberOfClumps()
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//             PUBLIC MEMBER FUNCTIONS

const forceModel& cfdemCloudMS::forceM(int i)
{
    return forceModel_[i];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
