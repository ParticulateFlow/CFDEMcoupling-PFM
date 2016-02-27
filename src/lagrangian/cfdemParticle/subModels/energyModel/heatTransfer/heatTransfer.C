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
#include "heatTransfer.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(heatTransfer, 0);

addToRunTimeSelectionTable(energyModel, heatTransfer, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
heatTransfer::heatTransfer
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    energyModel(dict,sm),
    LaEuScalarTemp(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    propsLaEuScalarTempDict_(dict.subDict("LaEuScalarTempProps")),
    Cp_(readScalar(propsLaEuScalarTempDict_.lookup("Cp")))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

heatTransfer::~heatTransfer()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void heatTransfer::energyContribution(volScalarField& EuField) const
{
    const volScalarField& rhoField = forceSubM(0).rhoField();
    temperatureContribution(EuField);
    EuField.internalField() *= rhoField.internalField()*Cp_;
}

void heatTransfer::temperatureContribution(volScalarField& EuField) const
{
    LaEuScalarTemp::manipulateScalarField(EuField);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
