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

Description
    Class scales force on particle close to boundary. Between inner and outer coordinate,
    scaling is linear between 1 and specified value, outside outer coordinate, scaling is
    constant with specified value.

Class
    scaleForce

SourceFiles
    scaleForce.C

\*---------------------------------------------------------------------------*/

#ifndef scaleForceBoundary_H
#define scaleForceBoundary_H

#include "forceSubModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scaleForceBoundary Declaration
\*---------------------------------------------------------------------------*/

class scaleForceBoundary
:
    public forceSubModel
{
private:

    dictionary propsDict_;

    scalar coordinateInner_;

    scalar coordinateOuter_;

    scalar dist_;

    scalar outerVal_;

    int dim_;

    int orientation_;

public:

    //- Runtime type information
    TypeName("scaleForceBoundary");


    // Constructors

        //- Construct from components
        scaleForceBoundary
        (
            const dictionary& dict,
            cfdemCloud& sm,
            forceModel& fm
        );

    // Destructor

        ~scaleForceBoundary();


    // Member Functions

        void partToArray(label, vector&, const vector&, const vector& Ufluid=vector::zero, scalar Cd=scalar(0)) const;

        word myType() const {return typeName; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
