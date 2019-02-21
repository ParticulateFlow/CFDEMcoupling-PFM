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
#include "weightSecondaryPhase.H"
#include "mathExtra.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(weightSecondaryPhase, 0);

addToRunTimeSelectionTable(otherForceModel,weightSecondaryPhase, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
weightSecondaryPhase::weightSecondaryPhase
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    otherForceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    volfracFieldName_(propsDict_.lookup("volfracFieldName")),
    alpha_(sm.mesh().lookupObject<volScalarField> (volfracFieldName_)),
    rho_("rho",dimensionSet(1,-3,0,0,0),0.0),
    g_("g",dimensionSet(0,1,-2,0,0),vector(0,0,-9.81))
{
    if (propsDict_.found("rho"))
        rho_.value()=readScalar(propsDict_.lookup ("rho"));
    else
        FatalError <<"Please specify density of secondary phase.\n" << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

weightSecondaryPhase::~weightSecondaryPhase()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volVectorField> weightSecondaryPhase::exportForceField()
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

    source = rho_ * alpha_ * g_;

    return tsource;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
