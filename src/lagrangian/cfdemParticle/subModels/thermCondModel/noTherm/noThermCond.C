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

#include "noThermCond.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noThermCond, 0);

addToRunTimeSelectionTable(thermCondModel, noThermCond, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
noThermCond::noThermCond
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    thermCondModel(dict,sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noThermCond::~noThermCond()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

tmp<volScalarField> noThermCond::thermCond() const
{
    tmp<volScalarField> tcond
    (
        new volScalarField
        (
            IOobject
            (
                "fake1",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
             ),
             particleCloud_.mesh(),
             dimensionedScalar
             (
              "zero",
              dimensionSet(1,1,-3,-1,0,0,0),
              0.0
             )
        )
    );

    return tcond;
}

tmp<volScalarField> noThermCond::thermDiff() const
{
    tmp<volScalarField> tdif
    (
        new volScalarField
        (
            IOobject
            (
                "fake2",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
             ),
             particleCloud_.mesh(),
             dimensionedScalar
             (
              "zero",
              dimensionSet(0,2,-1,0,0,0,0),
              0.0
             )
        )
    );

    return tdif;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
