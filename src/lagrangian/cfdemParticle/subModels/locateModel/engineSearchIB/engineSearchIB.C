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
#include "engineSearchIB.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(engineSearchIB, 0);

addToRunTimeSelectionTable
(
    engineSearch,
    engineSearchIB,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
engineSearchIB::engineSearchIB
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    engineSearch(dict.subDict(typeName + "Props"), sm),
    propsDict_(dict.subDict(typeName + "Props")),
    zSplit_(readLabel(propsDict_.lookup("zSplit"))),
    xySplit_(readLabel(propsDict_.lookup("xySplit"))),
    thetaSize_(180./zSplit_),
    phiSize_(360./xySplit_),
    deg2rad_(constant::mathematical::pi/180.),
    numberOfSatellitePoints_((zSplit_-1)*xySplit_+2)
{
    for (int countPoints = 0; countPoints < numberOfSatellitePoints_; ++countPoints)
    {
        satellitePoints_.push_back(generateSatellitePoint(countPoints));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

engineSearchIB::~engineSearchIB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


label engineSearchIB::findCell
(
    double** const& mask,
    double**& positions,
    int**& cellIDs,
    int size
) const
{
    bool checkPeriodicCells(particleCloud_.checkPeriodicCells());
    const boundBox& globalBb = particleCloud_.mesh().bounds();

    vector position;
    for (int index = 0; index < size; ++index)
    {
        cellIDs[index][0] = -1;
        double radius = particleCloud_.radius(index);

        if (radius > SMALL)
        {
            // create pos vector
            for (int i = 0; i < 3; i++) position[i] = positions[index][i];

            // find cell
            label oldID = cellIDs[index][0];
            cellIDs[index][0] = findSingleCell(position, oldID);

            if (cellIDs[index][0] < 0)
            {
                label altStartPos = -1;

                for (int countPoints = 0; countPoints < numberOfSatellitePoints_; ++countPoints)
                {
                    vector pos = getSatellitePoint(index, countPoints);

                    altStartPos = findSingleCell(pos,oldID);

                    //check for periodic domains
                    if (checkPeriodicCells)
                    {
                        for (int iDir = 0; iDir < 3; iDir++)
                        {
                            if (pos[iDir] > globalBb.max()[iDir])
                            {
                                pos[iDir] -= globalBb.max()[iDir] - globalBb.min()[iDir];
                            }
                            else if (pos[iDir] < globalBb.min()[iDir])
                            {
                                pos[iDir] += globalBb.max()[iDir] - globalBb.min()[iDir];
                            }
                        }

                        altStartPos = findSingleCell(pos, oldID);
                    }

                    if (altStartPos >= 0) // found position, we're done
                    {
                        cellIDs[index][0] = altStartPos;
                        break;
                    }
                }

            }
        }
    }

    return 1;
}

vector engineSearchIB::generateSatellitePoint(int countPoints) const
{
    // 1 point at bottom, 1 point at top
    if (countPoints == 0)
    {
        return vector(0., 0., 1.);
    }
    else if (countPoints == 1)
    {
        return vector(0., 0., -1.);
    }
    else
    {
        const scalar thetaLevel = (countPoints - 2) / xySplit_;
        const scalar theta = deg2rad_ * thetaSize_ * (thetaLevel + 1);
        const scalar phi = deg2rad_ * phiSize_ * (countPoints - 2 - thetaLevel * xySplit_);
        return vector(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
    }
}

vector engineSearchIB::getSatellitePoint(int index, int countPoints) const
{
    return particleCloud_.position(index)
         + particleCloud_.radius(index) * satellitePoints_[countPoints];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
