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
#include "gravity.H"
#include "mathExtra.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gravity, 0);

addToRunTimeSelectionTable(otherForceModel,gravity, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
gravity::gravity
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    otherForceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    rhoGFieldName_(propsDict_.lookupOrDefault<word>("rhoGFieldName","rho")),
    rhoG_(sm.mesh().lookupObject<volScalarField> (rhoGFieldName_)),
    g_("g",dimensionSet(0,1,-2,0,0),vector(0,0,-9.81))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gravity::~gravity()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volVectorField> gravity::exportForceField()
{
    tmp<volVectorField> tsource
    (
        new volVectorField
        (
            IOobject
            (
                "grav",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedVector
            (
                "zero",dimensionSet(1, -2, -2, 0, 0),vector::zero
            )
        )
    );

    volVectorField& source = tsource.ref();

    source = rhoG_ * voidfraction_ * g_;

    return tsource;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
