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

#include "scaleForceBoundary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(scaleForceBoundary, 0);

addToRunTimeSelectionTable
(
    forceSubModel,
    scaleForceBoundary,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
scaleForceBoundary::scaleForceBoundary
(
    const dictionary& dict,
    cfdemCloud& sm,
    forceModel& fm
)
:
    forceSubModel(dict,sm,fm),
    propsDict_(dict.subDict(typeName + "Props")),
    coordinateInner_(0.0),
    coordinateOuter_(0.0),
    outerVal_(0.0),
    orientation_(0)
{
    if (propsDict_.found("x1") && propsDict_.found("x2"))
    {
        coordinateInner_=readScalar(propsDict_.lookup ("x1"));
        coordinateOuter_=readScalar(propsDict_.lookup ("x2"));
        dim_ = 0;
        Info << "scaleForceBoundary: Limiting force in x direction."  << endl;
    }
    else if (propsDict_.found("y1") && propsDict_.found("y2"))
    {
        coordinateInner_=readScalar(propsDict_.lookup ("y1"));
        coordinateOuter_=readScalar(propsDict_.lookup ("y2"));
        dim_ = 1;
        Info << "scaleForceBoundary: Limiting force in y direction."  << endl;
    }
    else if (propsDict_.found("z1") && propsDict_.found("z2"))
    {
        coordinateInner_=readScalar(propsDict_.lookup ("z1"));
        coordinateOuter_=readScalar(propsDict_.lookup ("z2"));
        dim_ = 2;
        Info << "scaleForceBoundary: Limiting force in z direction."  << endl;
    }

    if (propsDict_.found("outerValue"))
        outerVal_=readScalar(propsDict_.lookup ("outerValue"));

    if(coordinateOuter_ > coordinateInner_)
        orientation_ = 1;
    else
        orientation_ = -1;

    dist_ = fabs(coordinateOuter_ - coordinateInner_);

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

scaleForceBoundary::~scaleForceBoundary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void scaleForceBoundary::partToArray
(
    label index,
    vector& dragTot,
    const vector& dragEx,
    const vector& Ufluid,
    scalar Cd
) const
{
    scalar scaleFac = 1.0;
    scalar weight = 1.0;
    scalar coordinate = particleCloud_.position(index).component(dim_);

    if (coordinate*orientation_ > coordinateInner_*orientation_)
    {
        weight = fabs( coordinateInner_ - coordinate) / dist_;
        if (weight > 1.0) weight = 1.0;
        scaleFac = (1.0 - weight) + weight * outerVal_;
    }

    dragTot *= scaleFac;
    Cd *= scaleFac;
    forceSubModel::partToArray(index,dragTot,dragEx,Ufluid,Cd);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
