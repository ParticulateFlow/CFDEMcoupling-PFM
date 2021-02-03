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

#include "terminalVelocity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(terminalVelocity, 0);

addToRunTimeSelectionTable
(
    forceModel,
    terminalVelocity,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
terminalVelocity::terminalVelocity
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
    mesh_(sm.mesh()),
    turbDissipationRateFieldName_(propsDict_.lookupOrDefault<word>("turbDissipationRateFieldName","")),
    turbDissipationRate_(NULL),
    existturbDissipationRateInObjReg_(false),
    liquidViscosity_(propsDict_.lookupOrDefault<scalar>("liquidViscosity", 1.0)),
    dragReductionFactor_(propsDict_.lookupOrDefault<scalar>("dragReductionFactor", 0.0)),
    terminalVel_(propsDict_.lookupOrDefault<scalar>("terminalVelocity", 0.0)),
    risingVelDir_(propsDict_.lookupOrDefault<scalar>("risingVelDir", 2))
{
    if(ignoreCellsName_ != "none")
    {
       ignoreCells_.set(new cellSet(particleCloud_.mesh(),ignoreCellsName_));
       Info << "isotropicFluctuations: ignoring rising velocity in cellSet " << ignoreCells_().name() <<
        " with " << ignoreCells_().size() << " cells." << endl;
    }
    else existIgnoreCells_ = false;

    if (turbDissipationRateFieldName_ != "")
    {
        existturbDissipationRateInObjReg_ = true;
        volScalarField& epsilon(const_cast<volScalarField&>(sm.mesh().lookupObject<volScalarField> (turbDissipationRateFieldName_)));
        turbDissipationRate_ = &epsilon;
    }
    else
    {
       turbDissipationRate_ = new volScalarField
        (
            IOobject
            (
                "turbDissipationRate",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0,2,-3,0,0), 0)
        );
    }
  
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

terminalVelocity::~terminalVelocity()
{
        if (!existturbDissipationRateInObjReg_) delete turbDissipationRate_;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

bool terminalVelocity::ignoreCell(label cell) const
{
    if (!existIgnoreCells_) return false;
    else return ignoreCells_()[cell];
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void terminalVelocity::setForce() const
{

    if (!existturbDissipationRateInObjReg_)
    {
        *turbDissipationRate_ = particleCloud_.turbulence().epsilon()();
    }
    vector position(0,0,0);
    scalar Urising=0.0;
    label cellI=-1;
    scalar radius=0.0;
    scalar epsilon=0.0;

    interpolationCellPoint<scalar> turbDissipationRateInterpolator_(*turbDissipationRate_);
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            Urising =0.0;
            if (cellI > -1 && !ignoreCell(cellI)) // particle found
            {
                if (interpolate_)
                {
                    position = particleCloud_.position(index);
                    epsilon = turbDissipationRateInterpolator_.interpolate(position,cellI);
                }
                else
                {
                    epsilon = (*turbDissipationRate_)[cellI];
                }

                    position = particleCloud_.position(index);
                    radius = particleCloud_.radius(index);

                    scalar dLambda = 2*radius*pow(epsilon,0.25)/pow(liquidViscosity_,0.75);

                    Urising = terminalVel_ / Foam::sqrt(1+(dragReductionFactor_*pow(dLambda,3)) );
                    particleCloud_.particleConvVels()[index][risingVelDir_] += Urising;
	    }
    }
   
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
