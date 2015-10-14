/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling
    
    Contributing authors:
    Thomas Lichtenegger, Gerhard Holzinger
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
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "forceModelRec.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(forceModelRec, 0);

defineRunTimeSelectionTable(forceModelRec, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
forceModelRec::forceModelRec
(
    const dictionary& dict,
    cfdemCloudRec& sm
)
:
    dict_(dict),
    particleCloud_(sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceModelRec::~forceModelRec()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //
void forceModelRec::partToArray
(
    label& index,
    const vector& dragTot
) const
{
        for(int j=0;j<3;j++) 
            DEMForces()[index][j] += dragTot[j];  
}

void forceModelRec::partToArrayU
(
    label& index,
    const vector& Ufluid
) const
{
        for(int j=0;j<3;j++)
           fluidVel()[index][j]=Ufluid[j]; 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
