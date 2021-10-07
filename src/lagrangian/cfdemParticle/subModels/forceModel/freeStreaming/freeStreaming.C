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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(freeStreaming, 0);

addToRunTimeSelectionTable
(
    forceModel,
    freeStreaming,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
freeStreaming::freeStreaming
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    interpolate_(propsDict_.lookupOrDefault<bool>("interpolation", false)),
    ignoreCellsName_(propsDict_.lookupOrDefault<word>("ignoreCellsName","none")),
    ignoreCells_(),
    existIgnoreCells_(true),
    UsRecFieldName_(propsDict_.lookupOrDefault<word>("granVelRecFieldName","UsRec")),
    UsRec_(sm.mesh().lookupObject<volVectorField> (UsRecFieldName_)),
    voidfractionRecFieldName_(propsDict_.lookupOrDefault<word>("voidfractionRecFieldName","voidfractionRec")),
    voidfractionRec_(sm.mesh().lookupObject<volScalarField> (voidfractionRecFieldName_)),
    critVoidfraction_(propsDict_.lookupOrDefault<scalar>("critVoidfraction", 1.0)),
    scalingFactor_(propsDict_.lookupOrDefault<scalar>("scalingFactor", 1.0)),
    particleDensity_(propsDict_.lookupOrDefault<scalar>("particleDensity",0.0)),
    gravAcc_(propsDict_.lookupOrDefault<vector>("g",vector(0.0,0.0,-9.81)))
{
    if(ignoreCellsName_ != "none")
    {
       ignoreCells_.set(new cellSet(particleCloud_.mesh(),ignoreCellsName_));
       Info << "isotropicFluctuations: ignoring fluctuations in cellSet " << ignoreCells_().name() <<
        " with " << ignoreCells_().size() << " cells." << endl;
    }
    else existIgnoreCells_ = false;
  
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

freeStreaming::~freeStreaming()
{}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

bool freeStreaming::ignoreCell(label cell) const
{
    if (!existIgnoreCells_) return false;
    else return ignoreCells_()[cell];
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void freeStreaming::setForce() const
{
    vector position(0,0,0);
    vector UNew(0,0,0);
    label cellI=0;
    scalar radius=0.0;
    scalar mass=0.0;
    vector grav(0,0,0);
    interpolationCellPoint<vector> UsRecInterpolator_(UsRec_);
    interpolationCellPoint<scalar> voidfractionRecInterpolator_(voidfractionRec_);
    
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            UNew =vector(0,0,0);
            if (cellI > -1 && !ignoreCell(cellI)) // particle found
            {
                // let particles in empty regions follow trajectories subject to gravity, otherwise stream
                if(voidfractionRec_[cellI] < critVoidfraction_)
                {
                    if( interpolate_ )
                    {
                        position = particleCloud_.position(index);
                        UNew = UsRecInterpolator_.interpolate(position,cellI);
                    }
                    else
                    {
                        UNew = UsRec_[cellI];
                    }
                    UNew /= scalingFactor_;
                }
                else
                {
                    position = particleCloud_.position(index);
                    radius = particleCloud_.radius(index);
                    mass = 4.188790205*radius*radius*radius * particleDensity_;
                    grav = mass*gravAcc_;
                    for(int j=0;j<3;j++)
                    {
                        particleCloud_.DEMForces()[index][j] += grav[j];
                    }
                    UNew=particleCloud_.velocity(index);
                }
                // write particle based data to global array
                for(int j=0;j<3;j++)
                {
                    particleCloud_.particleConvVels()[index][j] += UNew[j];
                }
	    }
    }
   
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
