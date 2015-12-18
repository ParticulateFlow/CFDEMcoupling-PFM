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
    fluctuationAmp_(readScalar(propsDict_.lookup("fluctuationAmp"))),
    ranGen_(osRandomInteger())
{
    // get the average cell side length to set reference fluctuation velocity
    scalar totalVolume(0.0);
    label numCells(0);
    scalar avLength(0.0);
    forAll(particleCloud_.mesh().cells(),cellI)
    {
        totalVolume += particleCloud_.mesh().V()[cellI];
	numCells++;
    }
    avLength=Foam::pow(totalVolume/numCells,1.0/3.0);
    refFluctuationVel_=avLength/particleCloud_.mesh().time().deltaTValue();
    Info << "\nUsing isotropic fluctuations with\n" << endl;
    Info << "\t\treference fluctuation velocity " << refFluctuationVel_ << " and\n" << endl;
    Info << "\t\tfluctuation amplification factor " << fluctuationAmp_ << ".\n" << endl;
}


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
    scalar pFluc(0.0);
   
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
		    pFluc=deltaVoidfrac/(1-voidfraction+SMALL);
		    if(deltaVoidfrac>0 && ranGen_.scalar01()<pFluc)
		    {
                        flucU=unitRndVec()*scaleFluctuations(deltaVoidfrac);
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
        fluctuation=refFluctuationVel_*fluctuationAmp_;
        return fluctuation;
    }
}

vector isotropicFluctuations::unitRndVec() const
{
    // algorithm according to:
    // Marsaglia. "Choosing a point from the surface of a sphere." The Annals of Mathematical Statistics 43.2 (1972): 645-646.
    scalar v1(0.0);
    scalar v2(0.0);
    scalar s(10.0);
    scalar s2(0.0);
    vector rvec(0,0,0);
    while(s>1.0)
    {
        v1=2*(ranGen_.scalar01()-0.5);
	v2=2*(ranGen_.scalar01()-0.5);
	s=v1*v1+v2*v2;
    }
    s2=Foam::sqrt(1-s);
    rvec[1]=2*v1*s2;
    rvec[2]=2*v2*s2;
    rvec[3]=1-2*s;
    return rvec;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
