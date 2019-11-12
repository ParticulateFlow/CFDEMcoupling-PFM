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
    kf0_
    (
        dimensionedScalar(transportProperties_.lookup("kf")).value()
    ),
    Cp_
    (
        dimensionedScalar(transportProperties_.lookup("Cp")).value()
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

energyModel::~energyModel()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

scalar energyModel::Cp() const
{
    return Cp_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //