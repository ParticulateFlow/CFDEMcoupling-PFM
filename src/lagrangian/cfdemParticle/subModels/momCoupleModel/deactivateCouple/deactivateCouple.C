/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling

    Contributing authors:
    Thomas Lichtenegger
    Copyright (C) 2015- Johannes Kepler University, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling academic.

    CFDEMcoupling academic is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling academic is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling academic.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "deactivateCouple.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(deactivateCouple, 0);

addToRunTimeSelectionTable
(
    momCoupleModel,
    deactivateCouple,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
deactivateCouple::deactivateCouple
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    momCoupleModel(dict,sm),
    fEmpty_
    (   IOobject
        (
            "fEmpty",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0), vector::zero),
        "zeroGradient"
    ),
    KslEmpty_
    (
        IOobject
        (
            "KslEmpty",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0), 0),
        "zeroGradient"
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

deactivateCouple::~deactivateCouple()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void deactivateCouple::resetMomSourceField()
{}

tmp<volVectorField> deactivateCouple::expMomSource()
{
    tmp<volVectorField> tsource(fEmpty_);
    return tsource;
}

tmp<volScalarField> deactivateCouple::impMomSource()
{
    tmp<volScalarField> tsource(KslEmpty_);
    return tsource;

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
