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

#include "ScaleForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ScaleForce, 0);

addToRunTimeSelectionTable
(
    forceSubModel,
    ScaleForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ScaleForce::ScaleForce
(
    const dictionary& dict,
    cfdemCloud& sm,
    forceModel& fm
)
:
    forceSubModel(dict,sm,fm),
    propsDict_(dict.subDict(typeName + "Props")),
    scaleFieldName_(propsDict_.lookup("ScaleField")),
    scaleField_(sm.mesh().lookupObject<volScalarField> (scaleFieldName_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ScaleForce::~ScaleForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void ScaleForce::partToArray
(
    label index,
    vector& dragTot,
    const vector& dragEx,
    const vector& Ufluid,
    scalar Cd
) const
{
    label cellI = particleCloud_.particleCell(index);
    dragTot *= scaleField_[cellI];
    forceSubModel::partToArray(index,dragTot,dragEx,Ufluid,Cd);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
