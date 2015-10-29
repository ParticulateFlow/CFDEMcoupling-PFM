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

#include "isotropicFluctuations.H"
#include "recModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(isotropicFluctuations, 0);

addToRunTimeSelectionTable
(
    forceModelRec,
    isotropicFluctuations,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
isotropicFluctuations::isotropicFluctuations
(
    const dictionary& dict,
    cfdemCloudRec& sm
)
:
    forceModelRec(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    interpolate_(propsDict_.lookupOrDefault<bool>("interpolation", false)),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    voidfractionRecFieldName_(propsDict_.lookupOrDefault<word>("voidfractionRecFieldName","voidfractionRec")),
    voidfractionRec_(sm.mesh().lookupObject<volScalarField> (voidfractionRecFieldName_)),
    critVoidfraction_(propsDict_.lookupOrDefault<scalar>("critVoidfraction", 1.0)),
    maxFluctuation_(readScalar(propsDict_.lookup("maxFluctuation"))),
    ranGen_(osRandomInteger())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

isotropicFluctuations::~isotropicFluctuations()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void isotropicFluctuations::setForce() const
{
    vector position(0,0,0);
    scalar voidfraction(0.0);
    scalar voidfractionRec(0.0);
    scalar deltaVoidfrac(0.0);

    vector flucU(0,0,0);
    label cellI=0;
   
    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<scalar> voidfractionRecInterpolator_(voidfractionRec_);
    
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
	    flucU =vector(0,0,0);
	    voidfraction=0.0;
	    voidfractionRec=0.0;
	    deltaVoidfrac=0.0;
            if (cellI > -1) // particle Found
            {
	        // particles in empty regions follow trajectories subject to gravity
		if(voidfractionRec_[cellI] < critVoidfraction_)
		{
                    if( interpolate_ )
                    {
                      position = particleCloud_.position(index);
		      voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
		      voidfractionRec = voidfractionRecInterpolator_.interpolate(position,cellI);
                    }
                    else
                    {
			voidfraction = voidfraction_[cellI];
			voidfractionRec = voidfractionRec_[cellI];
                    }
                    // write particle based data to global array
                    
                    deltaVoidfrac=voidfractionRec-voidfraction;
		    if(deltaVoidfrac>0)
		    {
		        for(int j=0;j<3;j++)
                            flucU[j]=2*(ranGen_.scalar01()-0.5);
		        flucU*=scaleFluctuations(deltaVoidfrac);
			partToArrayU(index,flucU);
		    }              
		}
	    }
    }
}

scalar isotropicFluctuations::scaleFluctuations(const scalar deltaVoidfrac) const
{
    scalar fluctuation;
    if(deltaVoidfrac<0.0)
    {
        return 0.0;
    }
    else
    {
        fluctuation=deltaVoidfrac*maxFluctuation_;
        return fluctuation;
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
