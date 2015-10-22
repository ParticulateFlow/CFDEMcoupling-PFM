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

#include "gradConcentration.H"
#include "recModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradConcentration, 0);

addToRunTimeSelectionTable
(
    forceModelRec,
    gradConcentration,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
gradConcentration::gradConcentration
(
    const dictionary& dict,
    cfdemCloudRec& sm
)
:
    forceModelRec(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    interpolate_(propsDict_.lookupOrDefault<bool>("interpolation", false)),
    scalingFactor_(propsDict_.lookupOrDefault<scalar>("scalingFactor", 1.0)),
    voidfractionRecFieldName_(propsDict_.lookupOrDefault<word>("voidfractionRecFieldName","voidfractionRec")),
    voidfractionRec_(sm.mesh().lookupObject<volScalarField> (voidfractionRecFieldName_)),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    deltaVoidfraction_
    (   IOobject
        (
            "deltaVoidfraction",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        voidfractionRec_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gradConcentration::~gradConcentration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gradConcentration::setForce() const
{
    deltaVoidfraction_ = voidfraction_-voidfractionRec_;
    volVectorField gradDeltaVoidfraction_ = fvc::grad(deltaVoidfraction_);
    vector position(0,0,0);
    vector force(0,0,0);
    label cellI=0;
    scalar scalDConc=0.0;
    interpolationCellPoint<scalar> dConcInterpolator_(deltaVoidfraction_);
    interpolationCellPoint<vector> gradConcInterpolator_(gradDeltaVoidfraction_);

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            force =vector(0,0,0);
            if (cellI > -1) // particle Found
            {

                if( interpolate_ )
                {
                  position = particleCloud_.position(index);
		  scalDConc=dConcInterpolator_.interpolate(position,cellI);
		  if(scalDConc<0)
                      force = gradConcInterpolator_.interpolate(position,cellI);
                }
                else
                {
		    if(deltaVoidfraction_[cellI]<0)
                        force = gradDeltaVoidfraction_[cellI];
                }
                // need a prefactor to relate concentration gradient to force
                 force*= scalingFactor_;
                
                // write particle based data to global array
                partToArray(index,force);
	    }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
