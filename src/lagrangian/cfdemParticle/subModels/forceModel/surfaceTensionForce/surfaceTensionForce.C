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

#include "surfaceTensionForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(surfaceTensionForce, 0);

addToRunTimeSelectionTable
(
    forceModel,
    surfaceTensionForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
surfaceTensionForce::surfaceTensionForce
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    stf_(sm.mesh().lookupObject<surfaceScalarField> ("surfaceTensionForce"))
{


    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();

    Info << "check if interpolation really works - use directly interpolationCellPoint<vector> ???" << endl;
    particleCloud_.checkCG(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

surfaceTensionForce::~surfaceTensionForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void surfaceTensionForce::setForce() const
{
    volVectorField reconstructedStf = fvc::reconstruct(stf_*particleCloud_.mesh().magSf());
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        //if(mask[index][0])
        //{
            // definition of spherical particle
            label cellI = particleCloud_.cellIDs()[index][0];
            scalar Vp = particleCloud_.particleVolume(index);
            if(cellI >-1.0) // particle found on proc domain
            {
                vector surfaceTensionForcep = Foam::vector(0,0,0);
                
                surfaceTensionForcep = Vp * reconstructedStf[cellI];

               // write particle based data to global array
               forceSubM(0).partToArray(index,surfaceTensionForcep,vector::zero);

            } // end if particle found on proc domain
        //}// end if in mask
    }// end loop particles
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
