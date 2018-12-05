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

#include "engineSearchMany2Many.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(engineSearchMany2Many, 0);

addToRunTimeSelectionTable
(
    locateModel,
    engineSearchMany2Many,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
engineSearchMany2Many::engineSearchMany2Many
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    engineSearch(dict.subDict(typeName + "Props"),sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

engineSearchMany2Many::~engineSearchMany2Many()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label engineSearchMany2Many::intersection
(
    const point& pStart,
    const point& pEnd
) const
{
    // find intersection with boundary
    label face = searchEngine_.findNearestBoundaryFace(pEnd);

    // try alternative
    if (face==-1)
    {
        face = searchEngine_.intersection(pStart,pEnd).index();

        if (face==-1 && mag(pStart-point(0,0,0))<SMALL)
        {
            point pStart2 = pEnd+0.0001*(pStart-pEnd)/mag(pStart-pEnd);
            face = searchEngine_.intersection(pStart2,pEnd).index();
        }
    }
    return face;
}

label engineSearchMany2Many::findCell
(
    double** const& mask,
    double**& positions,
    int**& cellIDs,
    int size
) const
{
    return 1; // locate is provided by many2may / one2one
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
