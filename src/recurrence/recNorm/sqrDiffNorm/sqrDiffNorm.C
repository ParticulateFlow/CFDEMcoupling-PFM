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
#include "sqrDiffNorm.H"
#include "recModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sqrDiffNorm, 0);

addToRunTimeSelectionTable
(
    recNorm,
    sqrDiffNorm,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
sqrDiffNorm::sqrDiffNorm
(
    const dictionary& dict,
    recBase& base
):
    diffNorm(dict, base, typeName)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sqrDiffNorm::~sqrDiffNorm()
{}

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * //



scalar sqrDiffNorm::normVSF(label ti, label tj)
{
    const volScalarField& t1( base_.recM().exportVolScalarField(fieldName_,ti) );
    const volScalarField& t2( base_.recM().exportVolScalarField(fieldName_,tj) );
    dimensionedScalar tNorm( fvc::domainIntegrate( sqr( t1 - t2 ) ) );

    return tNorm.value();
}

scalar sqrDiffNorm::normVVF(label ti, label tj)
{
    const volVectorField& t1( base_.recM().exportVolVectorField(fieldName_,ti) );
    const volVectorField& t2( base_.recM().exportVolVectorField(fieldName_,tj) );
    dimensionedScalar tNorm( fvc::domainIntegrate( magSqr( t1 - t2 ) ) );

    return tNorm.value();
}

scalar sqrDiffNorm::normSSF(label ti, label tj)
{
    const surfaceScalarField& t1( base_.recM().exportSurfaceScalarField(fieldName_,ti) );
    const surfaceScalarField& t2( base_.recM().exportSurfaceScalarField(fieldName_,tj) );
    volVectorField t12 (fvc::reconstruct( t1-t2 ) );
    dimensionedScalar tNorm( fvc::domainIntegrate( magSqr( t12 ) ) );

    return tNorm.value();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
