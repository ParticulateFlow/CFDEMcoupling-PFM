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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "energyModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(energyModel, 0);

defineRunTimeSelectionTable(energyModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
energyModel::energyModel
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            sm.mesh().time().constant(),
            sm.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    kf0Field_
    (
        IOobject
        (
            "kf0",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,1,-3,-1,0,0,0), 0.)
    ),
    CpField_
    (
        IOobject
        (
            "Cp",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,2,-2,-1,0,0,0), 0.)
    )
{
    // build constant fields for single phase case
    if (!particleCloud_.multiphase())
    {
        kf0Field_ = volScalarField
        (
            IOobject
            (
                "kf0",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedScalar(transportProperties_.lookup("kf"))
        );

        CpField_ = volScalarField
        (
            IOobject
            (
                "Cp",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedScalar(transportProperties_.lookup("Cp"))
        );
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

energyModel::~energyModel()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

const volScalarField& energyModel::kf0Field() const
{
    if (particleCloud_.multiphase())
    {
        return particleCloud_.mesh().lookupObject<volScalarField>("kf");
    }
    else
    {
        return kf0Field_;
    }
}

const volScalarField& energyModel::CpField() const
{
    if (particleCloud_.multiphase())
    {
        return particleCloud_.mesh().lookupObject<volScalarField>("Cp");
    }
    else
    {
        return CpField_;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //