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
//#include "recModel.H"
#include "forceModel.H"
namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template <class baseCloud>
cfdemCloudRec<baseCloud>::cfdemCloudRec
(
    const fvMesh& mesh
)
:
    baseCloud(mesh),
//     recModel_
//     (
//         recModel::New
//         (
//             couplingProperties_,
//             *this
//         )
//     ),
    coupleRecForce_(true)
{}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //
template <class baseCloud>
cfdemCloudRec<baseCloud>::~cfdemCloudRec()
{}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// shouldn't be needed anymore
// getDEMdata of base class sufficient
/*
void cfdemCloudRec::getDEMdata()
{
    dataExchangeM().getData("radius","scalar-atom",radii_);
    dataExchangeM().getData("x","vector-atom",positions_);
    dataExchangeM().getData("v","vector-atom",velocities_);
}
*/
template <class baseCloud>
void cfdemCloudRec<baseCloud>::giveDEMdata()
{
    // ab use fluidVel_ as difference between recurrence velocity and mean particle velocity
    baseCloud::dataExchangeM().giveData("vrec","vector-atom",baseCloud::fluidVel_);
    if(coupleRecForce_)
        baseCloud::dataExchangeM().giveData("dragforce","vector-atom",baseCloud::DEMForces_);
    if(baseCloud::verbose_) Info << "giveDEMdata done." << endl;
}

template <class baseCloud>
void cfdemCloudRec<baseCloud>::setForces()
{
    baseCloud::resetArray(baseCloud::fluidVel_,baseCloud::numberOfParticles(),3);
    baseCloud::resetArray(baseCloud::DEMForces_,baseCloud::numberOfParticles(),3);
    for (int i=0;i<baseCloud::nrForceModelsRec();i++) baseCloud::forceM(i).setForce();
}

template <class baseCloud>
void cfdemCloudRec<baseCloud>::setParticleForceField()
{}

// * * *   write top level fields   * * * //

// * * * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //
// void cfdemCloudRec::updateRecFields()
// {
//     recModel_->updateRecFields(); 
// }



// shouldn't be needed anymore
// evolve of base class sufficient, override setForces() with rec. forces
/*
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

        clockM().start(25,"dumpDEMdata");
        // do particle IO
        IOM().dumpDEMdata();
        clockM().stop("dumpDEMdata");
    }
    return doCouple;
}

*/



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
