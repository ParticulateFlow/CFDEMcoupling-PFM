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
    granKineticEnergy.C
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "granKineticEnergy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(granKineticEnergy, 0);

addToRunTimeSelectionTable
(
    forceModel,
    granKineticEnergy,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
granKineticEnergy::granKineticEnergy
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    vfluc_(NULL),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    granKineticEnergy_
    (   IOobject
        (
            "granKineticEnergy",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,2,-2,0,0), 0),
        "zeroGradient"
    )
{
    allocateMyArrays();
    granKineticEnergy_.write();


    // init force sub model
    setForceSubModels(propsDict_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

granKineticEnergy::~granKineticEnergy()
{
    particleCloud_.dataExchangeM().destroy(vfluc_,1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void granKineticEnergy::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal = 0.0;
    particleCloud_.dataExchangeM().allocateArray(vfluc_,initVal,1);
}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void granKineticEnergy::setForce() const
{
    allocateMyArrays();

    label cellI = 0;
    vector velfluc(0,0,0);


    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if(cellI >= 0)
        {
            velfluc = particleCloud_.velocity(index) - UsField_[cellI];
            vfluc_[index][0] = magSqr(velfluc);
        }
    }

    granKineticEnergy_.primitiveFieldRef() = 0.0;

    particleCloud_.averagingM().resetWeightFields();
    particleCloud_.averagingM().setScalarAverage
    (
        granKineticEnergy_,
        vfluc_,
        particleCloud_.particleWeights(),
        particleCloud_.averagingM().UsWeightField(),
        NULL
    );

    granKineticEnergy_ *= 0.5;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
