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

#include "turboEngineSearch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turboEngineSearch, 0);

addToRunTimeSelectionTable
(
    locateModel,
    turboEngineSearch,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
turboEngineSearch::turboEngineSearch
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    locateModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    treeSearch_(propsDict_.lookup("treeSearch")),
    bb_(particleCloud_.mesh().points(),false),
    searchEngine_(particleCloud_.mesh(), polyMesh::FACE_PLANES)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

turboEngineSearch::~turboEngineSearch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label turboEngineSearch::findCell
(
    double** const& mask,
    double**& positions,
    int**& cellIDs,
    int size
) const
{
    bool first=true;
    vector position;
    for(int index = 0;index < size; ++index)
    {
        //if(mask[index][0] && particleCloud_.radius(index) > SMALL)
        if(particleCloud_.radius(index) > SMALL)
        {
            // create pos vector
            for(int i=0;i<3;i++) position[i] = positions[index][i];

            // find cell
            if(first)
            {
                cellIDs[index][0] = searchEngine_.findCell(position,cellIDs[index][0],treeSearch_);
                first=false;
            }
            else
            {
                if(bb_.contains(position))
                    cellIDs[index][0] = searchEngine_.findCell(position,cellIDs[index][0],treeSearch_);
                else
                    cellIDs[index][0] = -1;
            }
        }
        else cellIDs[index][0] = -1;
    }

    return 1;
}

label turboEngineSearch::findSingleCell
(
    const vector& position,
    label oldCellID
) const
{
    // find cell
    return searchEngine_.findCell(position,oldCellID,treeSearch_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
