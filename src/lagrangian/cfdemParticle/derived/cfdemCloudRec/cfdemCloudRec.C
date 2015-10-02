/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling
    
    Contributing authors:
    Thomas Lichtenegger
    Copyright (C) 2015- Johannes Kepler University, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling academic.

    CFDEMcoupling academic is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling academic is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling academic.  If not, see <http://www.gnu.org/licenses/>.

Description
    This code is designed to realize coupled rCFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "cfdemCloudRec.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
cfdemCloudRec::cfdemCloudRec
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    recModel_
    (
        recModel::New
        (
            couplingProperties_,
            *this
        )
    )
    {}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //
cfdemCloudRec::~cfdemCloudRec()
{}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void cfdemCloudRec::getDEMdata()
{
    dataExchangeM().getData("radius","scalar-atom",radii_);
    dataExchangeM().getData("x","vector-atom",positions_);
    dataExchangeM().getData("v","vector-atom",velocities_);
}

void cfdemCloudRec::giveDEMdata()
{
    dataExchangeM().giveData("v","vector-atom",velocities_);
    dataExchangeM().giveData("dragforce","vector-atom",DEMForces_);
    if(verbose_) Info << "giveDEMdata done." << endl;
}

void cfdemCloud::setForces()
{
  //  resetArray(fluidVel_,numberOfParticles(),3);
  //  resetArray(impForces_,numberOfParticles(),3);
  //  resetArray(expForces_,numberOfParticles(),3);
  //  resetArray(DEMForces_,numberOfParticles(),3);
  //  resetArray(Cds_,numberOfParticles(),1);
  //  for (int i=0;i<cfdemCloud::nrForceModels();i++) cfdemCloud::forceM(i).setForce();
}

// * * *   write top level fields   * * * //

// * * * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * * //
void cfdemCloud::setForces()
{
  //  resetArray(fluidVel_,numberOfParticles(),3);
  //  resetArray(impForces_,numberOfParticles(),3);
  //  resetArray(expForces_,numberOfParticles(),3);
  //  resetArray(DEMForces_,numberOfParticles(),3);
  //  resetArray(Cds_,numberOfParticles(),1);
  //  for (int i=0;i<cfdemCloud::nrForceModels();i++) cfdemCloud::forceM(i).setForce();
}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void cfdemCloudRec::initRecFields()
{
    recM().initRecFields(); 
}

void cfdemCloudRec::updateRecFields()
{
    recM().updateRecFields(); 
}

bool cfdemCloudRec::evolve
(
    volScalarField& alpha,
    volVectorField& Us,
    volVectorField& U
)
{
    numberOfParticlesChanged_ = false;
    arraysReallocated_=false;
    bool doCouple=false;

    if(!ignore())
    {
        if (dataExchangeM().doCoupleNow())
        {
            Info << "\n Coupling..." << endl;
            dataExchangeM().couple(0);
            doCouple=true;

            // reset vol Fields
            clockM().start(16,"resetVolFields");
            if(verbose_)
            {
                Info << "couplingStep:" << dataExchangeM().couplingStep() 
                     << "\n- resetVolFields()" << endl;
            }
            averagingM().resetVectorAverage(averagingM().UsPrev(),averagingM().UsNext(),false);
            voidFractionM().resetVoidFractions();
            averagingM().resetWeightFields();
 
            if(verbose_) Info << "resetVolFields done." << endl;
            clockM().stop("resetVolFields");

            if(verbose_) Info << "- getDEMdata()" << endl;
            clockM().start(17,"getDEMdata");
            getDEMdata();
            clockM().stop("getDEMdata");
            if(verbose_) Info << "- getDEMdata done." << endl;

            // search cellID of particles
            clockM().start(18,"findCell");
            if(verbose_) Info << "- findCell()" << endl;
            findCells();
            if(verbose_) Info << "findCell done." << endl;
            clockM().stop("findCell");

            // set void fraction field
            clockM().start(19,"setvoidFraction");
            if(verbose_) Info << "- setvoidFraction()" << endl;
            voidFractionM().setvoidFraction(NULL,voidfractions_,particleWeights_,particleVolumes_,particleV_);
            if(verbose_) Info << "setvoidFraction done." << endl;
            clockM().stop("setvoidFraction");

            // set average particles velocity field
            clockM().start(20,"setVectorAverage");
            setVectorAverages();


            //Smoothen "next" fields            
            smoothingM().dSmoothing();
            smoothingM().smoothen(voidFractionM().voidFractionNext());
            
            clockM().stop("setVectorAverage");
        }
        
        //============================================
        //CHECK JUST TIME-INTERPOATE ALREADY SMOOTHENED VOIDFRACTIONNEXT AND UsNEXT FIELD 
        //      IMPLICIT FORCE CONTRIBUTION AND SOLVER USE EXACTLY THE SAME AVERAGED
        //      QUANTITIES AT THE GRID!
        Info << "\n timeStepFraction() = " << dataExchangeM().timeStepFraction() << endl;
        clockM().start(24,"interpolateEulerFields");

        // update voidFractionField
        alpha = voidFractionM().voidFractionInterp();
        if(dataExchangeM().couplingStep() < 2)
        {
            alpha.oldTime() = alpha; // supress volume src
            alpha.oldTime().correctBoundaryConditions();
        }
        alpha.correctBoundaryConditions();


        // update mean particle velocity Field
        Us = averagingM().UsInterp();
        Us.correctBoundaryConditions();

        clockM().stop("interpolateEulerFields");
        //============================================

        if(doCouple)
        {
            // set particles forces
            clockM().start(21,"setForce");
            if(verbose_) Info << "- setForce(forces_)" << endl;
            setForces();
            if(verbose_) Info << "setForce done." << endl;
            clockM().stop("setForce");

            // get particles velocities
            clockM().start(22,"setParticlesVelocities");
            if(verbose_) Info << "- setParticlesVelocities" << endl;
   //         setParticleForceField();
            if(verbose_) Info << "- setParticlesVelocities done." << endl;
            clockM().stop("setParticlesVelocities");

            // write DEM data
            if(verbose_) Info << " -giveDEMdata()" << endl;
            clockM().start(23,"giveDEMdata");
            giveDEMdata();
            clockM().stop("giveDEMdata");

            dataExchangeM().couple(1);
        }


        if(verbose_){
            #include "debugInfo.H"
        }

        clockM().start(25,"dumpDEMdata");
        // do particle IO
        IOM().dumpDEMdata();
        clockM().stop("dumpDEMdata");
    }
    return doCouple;
}





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
