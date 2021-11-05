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

#include "noDiffCoeff.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noDiffCoeff, 0);

addToRunTimeSelectionTable(diffCoeffModel, noDiffCoeff, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
noDiffCoeff::noDiffCoeff
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    diffCoeffModel(dict,sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noDiffCoeff::~noDiffCoeff()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

tmp<volScalarField> noDiffCoeff::diffCoeff() const
{
    tmp<volScalarField> dcoeff
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
              dimensionSet(0,2,-1,0,0,0,0),
              0.0
             )
        )
    );

    return dcoeff;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
