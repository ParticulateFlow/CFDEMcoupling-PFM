/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling

    Contributing authors:
    Thomas Lichtenegger
    Copyright (C) 2015- Johannes Kepler University, Linz

    Paul Kieckhefen, TUHH
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
#include "readNorm.H"
#include "recModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(readNorm, 0);

addToRunTimeSelectionTable
(
    recNorm,
    readNorm,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
readNorm::readNorm
(
    const dictionary& dict,
    recBase& base
)
:
    recNorm(dict, base),
    propsDict_(dict.subDict(typeName + "Props")),
    recMatName_(propsDict_.lookupOrDefault<word>("recMatName", "recurrenceMatrix"))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

readNorm::~readNorm()
{}

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * //

void readNorm::computeRecMatrix()
{
    Info << nl << type() << ": reading recurrence matrix " << recMatName_ << nl << endl;
    SymmetricSquareMatrix<scalar>& recurrenceMatrix( base_.recM().recurrenceMatrix() );
    if (Pstream::master())
    {
        IFstream matrixFile(recMatName_);
        matrixFile >> recurrenceMatrix;
    }
    Pstream::scatter(recurrenceMatrix);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
