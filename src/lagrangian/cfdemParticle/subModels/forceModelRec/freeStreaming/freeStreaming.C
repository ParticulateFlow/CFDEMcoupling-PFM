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
    allowFluctuations_(propsDict_.lookupOrDefault<bool>("fluctuations", false)),
    UsFieldName_(propsDict_.lookupOrDefault<word>("granVelFieldName","Us")),
    Us_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    UsRecFieldName_(propsDict_.lookupOrDefault<word>("granVelRecFieldName","UsRec")),
    UsRec_(sm.mesh().lookupObject<volVectorField> (UsRecFieldName_)),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    voidfractionRecFieldName_(propsDict_.lookupOrDefault<word>("voidfractionRecFieldName","voidfractionRec")),
    voidfractionRec_(sm.mesh().lookupObject<volScalarField> (voidfractionRecFieldName_)),
    critVoidfraction_(propsDict_.lookupOrDefault<scalar>("critVoidfraction", 1.0)),
    particleDensity_(propsDict_.lookupOrDefault<scalar>("particleDensity",0.0)),
    gravAcc_(propsDict_.lookupOrDefault<vector>("g",vector(0.0,0.0,-9.81)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

freeStreaming::~freeStreaming()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void freeStreaming::setForce() const
{
    vector position(0,0,0);
    vector Us(0,0,0);
    vector UsRec(0,0,0);
    scalar voidfraction(0.0);
    scalar voidfractionRec(0.0);
    vector newUs(0,0,0);
    vector flucU(0,0,0);
    label cellI=0;
    scalar radius=0.0;
    scalar mass=0.0;
    vector grav(0,0,0);
    interpolationCellPoint<vector> UsInterpolator_(Us_);
    interpolationCellPoint<vector> UsRecInterpolator_(UsRec_);
    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<scalar> voidfractionRecInterpolator_(voidfractionRec_);
    
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            Us =vector(0,0,0);
	    UsRec =vector(0,0,0);
	    voidfraction=0.0;
	    voidfractionRec=0.0;
            if (cellI > -1) // particle Found
            {
	        // let particles in empty regions follow trajectories subject to gravity, ohterwise stream
		if(voidfractionRec_[cellI] < critVoidfraction_)
		{
                    if( interpolate_ )
                    {
                      position = particleCloud_.position(index);
                      Us = UsInterpolator_.interpolate(position,cellI);
		      UsRec = UsRecInterpolator_.interpolate(position,cellI);
		      voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
		      voidfractionRec = voidfractionRecInterpolator_.interpolate(position,cellI);
                    }
                    else
                    {
		        Us = Us_[cellI];
                        UsRec = UsRec_[cellI];
			voidfraction = voidfraction_[cellI];
			voidfractionRec = voidfractionRec_[cellI];
                    }
                    // write particle based data to global array
                    newUs=UsRec;
		    if(allowFluctuations_)
		    {
		        for(int j=0;j<3;j++)
                            flucU[j]=particleVel()[index][j]-Us;
		        newUs+=scaleFluctuations(voidfractionRec,voidfraction)*flucU;
		    }
                    partToArrayU(index,newUs);
		}
		else
		{
		    radius = particleCloud_.radius(index);
                    mass = 4.188790205*radius*radius*radius * particleDensity_;
		    grav = mass*gravAcc_;	    
		    partToArray(index,grav);
		    for(int j=0;j<3;j++)
                            newUs[j]=particleVel()[index][j]
		    partToArrayU(index,newUs);
		}
	    }
    }
}


scalar freeStreaming::scaleFluctuations(const scalar voidfracRec, const scalar voidfrac)
{
    scalar deltaVoidfrac=voidfracRec-voidfrac;
    if(deltaVoidfrac<0.0)
    {
        return 0.0;
    }
    else
    {
        return deltaVoidfrac;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
