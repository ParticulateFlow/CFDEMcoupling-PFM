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

#include "dividedVoidFraction.H"
#include "addToRunTimeSelectionTable.H"
#include "locateModel.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dividedVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    dividedVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dividedVoidFraction::dividedVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(propsDict_.found("verbose")),
    procBoundaryCorrection_(propsDict_.lookupOrDefault<Switch>("procBoundaryCorrection", false)),
    alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
    alphaLimited_(0),
    tooMuch_(0.0),
    interpolation_(propsDict_.found("interpolation")),
    cfdemUseOnly_(propsDict_.lookupOrDefault<bool>("cfdemUseOnly", false))
{
    maxCellsPerParticle_ = numberOfMarkerPoints;

    if (alphaMin_ > 1.0 || alphaMin_ < 0.01)
        Warning << "alphaMin should be < 1 and > 0.01 !!!" << endl;

    if (interpolation_)
    {
        Warning << "interpolation for dividedVoidFraction does not yet work correctly!" << endl;
        Info << "Using interpolated voidfraction field - do not use this in combination with interpolation in drag model!" << endl;
    }

    checkWeightNporosity(propsDict_);

    if (procBoundaryCorrection_)
    {
        if (!(particleCloud_.locateM().type() == "engineIB"))
        {
            FatalError << typeName << ": You are requesting procBoundaryCorrection, this requires the use of engineIB!\n"
                       << abort(FatalError);
        }
    }
    else
    {
        if (particleCloud_.locateM().type() == "engineIB")
        {
            FatalError << typeName << ": You are using engineIB, this requires using procBoundaryCorrection=true!\n"
                       << abort(FatalError);
        }
    }

    // generate marker points on unit sphere
    label m = 0;
    offsets[m][0] = offsets[m][1] = offsets[m][2] = 0.0;
    ++m;

    // for 2 different radii
    scalar r1 = cbrt(1.0/numberOfMarkerPoints);
    scalar r2 = cbrt(15.0/numberOfMarkerPoints);
    scalar r[] = { 0.75 * (r2*r2*r2*r2 - r1*r1*r1*r1)/(r2*r2*r2 - r1*r1*r1),
                   0.75 * (1.0 - r2*r2*r2*r2)/(1.0 - r2*r2*r2) };

    for (label ir = 0; ir < 2; ++ir)
    {
        // try 8 subpoints derived from spherical coordinates
        for (scalar zeta = M_PI_4; zeta < constant::mathematical::twoPi; zeta += constant::mathematical::piByTwo)
        {
            for (scalar theta = M_PI_4; theta < constant::mathematical::pi; theta += constant::mathematical::piByTwo)
            {
                offsets[m][0] = r[ir] * Foam::sin(theta) * Foam::cos(zeta);
                offsets[m][1] = r[ir] * Foam::sin(theta) * Foam::sin(zeta);
                offsets[m][2] = r[ir] * Foam::cos(theta);
                ++m;
            }
        }
        // try 2 more subpoints for each coordinate direction (6 total)
        for (label j = -1; j <= 1; j += 2)
        {
            offsets[m][0] = r[ir] * j;
            offsets[m][1] = 0.;
            offsets[m][2] = 0.;
            ++m;

            offsets[m][0] = 0.;
            offsets[m][1] = r[ir] * j;
            offsets[m][2] = 0.;
            ++m;

            offsets[m][0] = 0.;
            offsets[m][1] = 0.;
            offsets[m][2] = r[ir] * j;
            ++m;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dividedVoidFraction::~dividedVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dividedVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes, double**& particleV)
{
    if (cfdemUseOnly_)
        reAllocArrays(particleCloud_.numberOfParticles());
    else
        reAllocArrays();

    vector position(0.,0.,0.);
    label cellID = -1;
    scalar radius(-1.);
    scalar volume(0.);
    scalar cellVol(0.);
    scalar scaleVol = weight();
    scalar scaleRadius = pow(porosity(),1./3.);
    const boundBox& globalBb = particleCloud_.mesh().bounds();

    for (int index=0; index < particleCloud_.numberOfParticles(); index++)
    {
        if (!checkParticleType(index)) continue; //skip this particle if not correct type

        //if(mask[index][0])
        //{
            // reset

            for (int subcell=0; subcell < cellsPerParticle_[index][0]; subcell++)
            {
                particleWeights[index][subcell] = 0.;
                particleVolumes[index][subcell] = 0.;
            }
            particleV[index][0] = 0.;

            cellsPerParticle_[index][0] = 1;
            position = particleCloud_.position(index);
            cellID = particleCloud_.cellIDs()[index][0];
            radius = particleCloud_.radius(index);
            if (multiWeights_) scaleVol = weight(index);
            volume = Vp(index,radius,scaleVol);
            radius *= scaleRadius;
            cellVol = 0;

            //--variables for sub-search
            int nPoints = numberOfMarkerPoints;
            int nNotFound = 0, nUnEqual = 0, nTotal = 0;
            vector offset(0.,0.,0.);
            int cellsSet = 0;

            if (procBoundaryCorrection_)
            {
                label cellWithCenter(-1);
                // switch off cellIDs for force calc if steming from parallel search success
                cellWithCenter = particleCloud_.locateM().findSingleCell(position,cellID);
                particleCloud_.cellIDs()[index][0] = cellWithCenter;
            }

            if (cellID >= 0)  // particel centre is in domain
            {
                cellVol = particleCloud_.mesh().V()[cellID];

                if (procBoundaryCorrection_)
                {
                    offset = radius * offsets[0];
                    #include "setWeightedSource.H"   // set source terms at position+offset
                }

                for (label i = 1; i < numberOfMarkerPoints; ++i)
                {
                    offset = radius * offsets[i];
                    #include "setWeightedSource.H"   // set source terms at position+offset
                }

                if (cellsSet > maxCellsPerParticle_ || cellsSet < 0)
                {
                    Info << "ERROR  cellsSet =" << cellsSet << endl;
                }

                if (!procBoundaryCorrection_)
                {
                    // set source for particle center; source 1/nPts+weight of all subpoints that have not been found
                    scalar centreWeight = 1./nPoints*(nPoints-cellsSet);

                    // update voidfraction for each particle read
                    scalar newAlpha = voidfractionNext_[cellID]- volume*centreWeight/cellVol;
                    if (newAlpha > alphaMin_)
                    {
                        voidfractionNext_[cellID] = newAlpha;
                    }
                    else
                    {
                        voidfractionNext_[cellID] = alphaMin_;
                        tooMuch_ += (alphaMin_-newAlpha) * cellVol;
                    }

                    // store cellweight for each particle --- this should be done for subpoints as well!!
                    particleWeights[index][0] += centreWeight;

                    // store particleVolume for each particle
                    particleVolumes[index][0] += volume*centreWeight;
                    particleV[index][0] += volume*centreWeight;
                }
            }// end if in cell
        //}// end if in mask
    }// end loop all particles
    voidfractionNext_.correctBoundaryConditions();

    // reset counter of lost volume
    if (verbose_) Pout << "Total particle volume neglected: " << tooMuch_ << endl;
    tooMuch_ = 0.;

    // bring voidfraction from Eulerian Field to particle array
    //interpolationCellPoint<scalar> voidfractionInterpolator_(voidfractionNext_);
    //scalar voidfractionAtPos(0);
    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {
        {
            for (int subcell=0; subcell < cellsPerParticle_[index][0]; subcell++)
            {
                label cellID = particleCloud_.cellIDs()[index][subcell];

                if (cellID >= 0)
                {
                    voidfractions[index][subcell] = voidfractionNext_[cellID];
                }
                else
                {
                    voidfractions[index][subcell] = -1.;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
