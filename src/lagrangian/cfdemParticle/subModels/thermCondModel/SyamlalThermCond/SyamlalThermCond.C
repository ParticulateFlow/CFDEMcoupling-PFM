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
#include "SyamlalThermCond.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SyamlalThermCond, 0);

addToRunTimeSelectionTable
(
    thermCondModel,
    SyamlalThermCond,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
SyamlalThermCond::SyamlalThermCond
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    thermCondModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SyamlalThermCond::~SyamlalThermCond()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void SyamlalThermCond::calcThermCond()
{
    forAll(thermCondField_,cellI)
    {
        if (1-voidfraction_[cellI] < SMALL) thermCondField_[cellI] = kf0_.value();
        else if (voidfraction_[cellI] < SMALL) thermCondField_[cellI] = 0.0;
        else thermCondField_[cellI] = (1-sqrt(1-voidfraction_[cellI]+SMALL)) / (voidfraction_[cellI]) * kf0_.value();
    }

    thermCondField_.correctBoundaryConditions();

    calcBoundaryCorrections();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
