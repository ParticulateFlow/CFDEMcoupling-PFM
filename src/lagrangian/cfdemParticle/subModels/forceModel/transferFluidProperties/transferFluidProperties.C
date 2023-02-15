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
    propsDict_(dict.subDict(typeName + "Props"))
{
    particleCloud_.registerParticleProperty<double**>("fluidDensity",1);
    particleCloud_.registerParticleProperty<double**>("fluidViscosity",1);

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_VERBOSE,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(SW_INTERPOLATION,true); // activate search for interpolate switch
    forceSubM(0).readSwitches();
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

    interpolationCellPoint<scalar> rhoInterpolator_(rhoField);
    interpolationCellPoint<scalar> nufInterpolator_(nufField);

    label cellI = 0;
    double rho = 0.;
    double nuf = 0.;
    vector position(0,0,0);

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if (cellI >= 0)
        {
            if(forceSubM(0).interpolation())
            {
                position = particleCloud_.position(index);
                rho = rhoInterpolator_.interpolate(position,cellI);
                nuf = nufInterpolator_.interpolate(position,cellI);
            }
            else
            {
                rho = rhoField[cellI];
                nuf = nufField[cellI];
            }

            fluidDensity_[index][0]   = rho;
            fluidViscosity_[index][0] = nuf*rho;
        }
    }

    particleCloud_.dataExchangeM().giveData("fluidDensity","scalar-atom",fluidDensity_);
    particleCloud_.dataExchangeM().giveData("fluidViscosity","scalar-atom",fluidViscosity_);
    
    if (forceSubM(0).verbose()) Info << "give data done" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
