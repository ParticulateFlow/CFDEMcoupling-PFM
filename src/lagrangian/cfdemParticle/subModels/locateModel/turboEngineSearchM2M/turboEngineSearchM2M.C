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

#include "turboEngineSearchM2M.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turboEngineSearchM2M, 0);

addToRunTimeSelectionTable
(
    locateModel,
    turboEngineSearchM2M,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
turboEngineSearchM2M::turboEngineSearchM2M
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    turboEngineSearch(dict.subDict(typeName + "Props"),sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

turboEngineSearchM2M::~turboEngineSearchM2M()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label turboEngineSearchM2M::findCell
(
    double** const& mask,
    double**& positions,
    double**& cellIDs,
    int size
) const
{
    // search should already be done by M2M

    return 1;
}

label turboEngineSearchM2M::intersection
(
    const point& pStart,
    const point& pEnd
) const
{
    // find intersection with boundary
    //pointIndexHit hit=searchEngine_.intersection(pStart,pEnd);
    //Info << "hit.index()=" << hit.index()<< endl;
    
    //return searchEngine_.findNearestBoundaryFace(pStart);
    return searchEngine_.intersection(pStart,pEnd).index();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
