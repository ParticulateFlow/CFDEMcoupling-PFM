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
    transfer fluid properties to LIGGGHTS

SourceFiles
    transferFluidProperties.C
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "transferFluidProperties.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(transferFluidProperties, 0);

addToRunTimeSelectionTable
(
    forceModel,
    transferFluidProperties,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
transferFluidProperties::transferFluidProperties
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false))
{
    particleCloud_.registerParticleProperty<double**>("fluidDensity",1);
    particleCloud_.registerParticleProperty<double**>("fluidViscosity",1);

    // init force sub model
    setForceSubModels(propsDict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

transferFluidProperties::~transferFluidProperties()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void transferFluidProperties::setForce() const
{
    double**& fluidDensity_   = particleCloud_.getParticlePropertyRef<double**>("fluidDensity");
    double**& fluidViscosity_ = particleCloud_.getParticlePropertyRef<double**>("fluidViscosity");

    const volScalarField& rhoField = forceSubM(0).rhoField();
    const volScalarField& nufField = forceSubM(0).nuField();

    label cellI = 0;

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if (cellI >= 0)
        {
            fluidDensity_[index][0] = rhoField[cellI];
            fluidViscosity_[index][0] = nufField[cellI] * rhoField[cellI];
        }
    }

    particleCloud_.dataExchangeM().giveData("fluidDensity","scalar-atom",fluidDensity_);
    particleCloud_.dataExchangeM().giveData("fluidViscosity","scalar-atom",fluidViscosity_);
    
    if (verbose_) Info << "give data done" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
