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
#include "secondaryPhaseInducedBuoyancy.H"
#include "mathExtra.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(secondaryPhaseInducedBuoyancy, 0);

addToRunTimeSelectionTable(otherForceModel, secondaryPhaseInducedBuoyancy, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
secondaryPhaseInducedBuoyancy::secondaryPhaseInducedBuoyancy
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    otherForceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    diffRho_("diffRho",dimensionSet(1,-3,0,0,0),0.0),
    g_("g",dimensionSet(0,1,-2,0,0),vector(0,0,-9.81))
{
    if (propsDict_.found("diffRho"))
        diffRho_.value()=readScalar(propsDict_.lookup ("diffRho"));
    else
        FatalError <<"Please specify density difference of the phases.\n" << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

secondaryPhaseInducedBuoyancy::~secondaryPhaseInducedBuoyancy()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volVectorField> secondaryPhaseInducedBuoyancy::exportForceField()
{
    tmp<volVectorField> tsource
    (
        new volVectorField
        (
            IOobject
            (
                "secondaryPhaseInducedBuoyancy",
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

    source = diffRho_ * (1.0-voidfraction_) * g_;

    return tsource;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
