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
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "freeStreaming.H"
#include "recModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(freeStreaming, 0);

addToRunTimeSelectionTable
(
    forceModelRec,
    freeStreaming,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
freeStreaming::freeStreaming
(
    const dictionary& dict,
    cfdemCloudRec& sm
)
:
    forceModelRec(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    interpolate_(propsDict_.lookupOrDefault<bool>("interpolation", false)),
    Usrec_(particleCloud_.recM().Us())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

freeStreaming::~freeStreaming()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void freeStreaming::setForce() const
{
    vector position(0,0,0);
    vector Us(0,0,0);
    label cellI=0;
    interpolationCellPoint<vector> UInterpolator_(Usrec_);
   
    // dummy variables
    vector drag(0,0,0);

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            Us =vector(0,0,0);
            if (cellI > -1) // particle Found
            {

                if( interpolate_ )
                {
                  position = particleCloud_.position(index);
                  Us = UInterpolator_.interpolate(position,cellI);
                }
                else
                {
                    Us = Usrec_[cellI];
                }
            // write particle based data to global array
                partToArray(index,drag,Us);
	    }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
