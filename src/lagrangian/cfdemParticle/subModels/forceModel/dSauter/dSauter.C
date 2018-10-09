/*---------------------------------------------------------------------------*\
License
    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2015- Thomas Lichtenegger, JKU Linz, Austria

Description
    calculates the Sauter mean diameter \sum d_i^3 / \sum d_i^2

SourceFiles
    dSauter.C
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "dSauter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dSauter, 0);

addToRunTimeSelectionTable
(
    forceModel,
    dSauter,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dSauter::dSauter
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    multiTypes_(false),
    d2_(NULL),
    d3_(NULL),
    typeCG_(propsDict_.lookupOrDefault<scalarList>("coarseGrainingFactors",scalarList(1,1.0))),
    d2Field_
    (   IOobject
        (
            "d2Field",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,2,0,0,0), 0)
    ),
    d3Field_
    (   IOobject
        (
            "d3Field",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,3,0,0,0), 0)
    ),
    dSauter_
    (   IOobject
        (
            "dSauter",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,1,0,0,0), 0),
        "zeroGradient"
    )
{
    if (typeCG_.size()>1) multiTypes_ = true;
    allocateMyArrays();
    dSauter_.write();


    // init force sub model
    setForceSubModels(propsDict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dSauter::~dSauter()
{
    particleCloud_.dataExchangeM().destroy(d2_,1);
    particleCloud_.dataExchangeM().destroy(d3_,1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

void dSauter::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal = 0.0;
    particleCloud_.dataExchangeM().allocateArray(d2_,initVal,1);  // field/initVal/with/lenghtFromLigghts
    particleCloud_.dataExchangeM().allocateArray(d3_,initVal,1);
}

// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void dSauter::setForce() const
{
    if (typeCG_.size()>1 || typeCG_[0] > 1.0)
    {
        Info << "dSauter using CG factor(s) = " << typeCG_ << endl;
    }

    allocateMyArrays();

    label cellI = 0;
    label partType = 1;
    scalar cg = typeCG_[0];
    scalar ds = 0.0;
    scalar effVolFac = 1.0;


    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if (cellI >= 0)
        {
            if (particleCloud_.getParticleEffVolFactors()) 
            {
                effVolFac = particleCloud_.particleEffVolFactor(index);
            }
            if (multiTypes_) 
            {
                partType = particleCloud_.particleType(index);
                cg = typeCG_[partType - 1];
            }
            ds = particleCloud_.d(index);
            d2_[index][0] = ds*ds*effVolFac*cg;
            d3_[index][0] = ds*ds*ds*effVolFac;
        }
    }

    d2Field_.primitiveFieldRef() = 0.0;
    d3Field_.primitiveFieldRef() = 0.0;

    particleCloud_.averagingM().setScalarSum
    (
        d2Field_,
        d2_,
        particleCloud_.particleWeights(),
        NULL
    );

    particleCloud_.averagingM().setScalarSum
    (
        d3Field_,
        d3_,
        particleCloud_.particleWeights(),
        NULL
    );

    forAll(dSauter_,cellI)
    {
        if (d2Field_[cellI] > ROOTVSMALL)
        {
            dSauter_[cellI] = d3Field_[cellI] / d2Field_[cellI];
        }
        else
        {
            dSauter_[cellI] = SMALL;
        }
    }

    dSauter_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
