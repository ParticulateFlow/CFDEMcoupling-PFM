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

#include "error.H"
#include "voidFractionModel.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(voidFractionModel, 0);

defineRunTimeSelectionTable(voidFractionModel, dictionary);

// * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
voidFractionModel::voidFractionModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    multiWeights_(false),
    voidfractionPrev_
    (   IOobject
        (
            "voidfractionPrev",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volScalarField> ("voidfraction")
        /*sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 1)*/
    ),
    voidfractionNext_
    (   IOobject
        (
            "voidfractionNext",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volScalarField> ("voidfraction")
        /*sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 1)*/
    ),
    cellsPerParticle_(NULL),
    maxCellsPerParticle_(1),
    weight_(1.),
    porosity_(1.)
{
    particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,1);
    if (particleCloud_.getParticleEffVolFactors()) multiWeights_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

voidFractionModel::~voidFractionModel()
{
    particleCloud_.dataExchangeM().destroy(cellsPerParticle_,1);
}

// * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //
tmp<volScalarField> voidFractionModel::voidFractionInterp() const
{
    const scalar tsf = particleCloud_.dataExchangeM().timeStepFraction();

    return tmp<volScalarField>
    (
        new volScalarField("alpha_voidFractionModel", (1. - tsf) * voidfractionPrev_ + tsf * voidfractionNext_)
    );
}

void voidFractionModel::resetVoidFractions()
{
    voidfractionPrev_.ref() = voidfractionNext_.ref();
    voidfractionNext_.ref() = 1.;
}

int** const& voidFractionModel::cellsPerParticle() const
{
    return cellsPerParticle_;
}

int voidFractionModel::maxCellsPerParticle() const
{
    return maxCellsPerParticle_;
}

void voidFractionModel::reAllocArrays()
{
    if(particleCloud_.numberOfParticlesChanged())
    {
        // get arrays of new length
        particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,1);
    }
}

void voidFractionModel::reAllocArrays(int nP)
{
    if(particleCloud_.numberOfParticlesChanged())
    {
        // get arrays of new length
        particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,1,nP);
    }
}

scalar voidFractionModel::pointInParticle(int index, const vector& positionCenter, const vector& point, double scale) const
{
    const scalar radius = particleCloud_.radius(index);

    if(radius > SMALL)
    {
        scalar pointDistSq = magSqr(point - positionCenter);
        return pointDistSq / (scale*scale * radius*radius) - 1.0;
    }
    else
    {
        return 0.;
    }
}

//Function to determine minimal distance of point
//to one of the periodic images of a particle
scalar voidFractionModel::minPeriodicDistance(int index,
                                           const vector& cellCentrePosition,
                                           const vector& positionCenter,
                                           const boundBox& globalBb,
                                           vector& minPeriodicPos) const
{
    scalar f = VGREAT;
    vector positionCenterPeriodic;

    for(label xDir=-1; xDir<=1; ++xDir)
    {
        positionCenterPeriodic[0] =  positionCenter[0]
                                  + static_cast<scalar>(xDir)
                                  * (globalBb.max()[0]-globalBb.min()[0]);
        for(label yDir=-1; yDir<=1; ++yDir)
        {
            positionCenterPeriodic[1] =  positionCenter[1]
                                      + static_cast<scalar>(yDir)
                                      * (globalBb.max()[1]-globalBb.min()[1]);
            for(label zDir=-1; zDir<=1; ++zDir)
            {
                positionCenterPeriodic[2] =  positionCenter[2]
                                          + static_cast<scalar>(zDir)
                                          * (globalBb.max()[2]-globalBb.min()[2]);

                if(pointInParticle(index, positionCenterPeriodic, cellCentrePosition) < f)
                {
                    f = pointInParticle(index, positionCenterPeriodic, cellCentrePosition);
                    minPeriodicPos = positionCenterPeriodic;
                }
            }
        }
    }

    return f;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
