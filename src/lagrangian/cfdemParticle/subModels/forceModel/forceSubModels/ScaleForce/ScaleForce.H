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


Class
    ScaleForce

SourceFiles
    ScaleForce.C

\*---------------------------------------------------------------------------*/

#ifndef ScaleForce_H
#define ScaleForce_H

#include "forceSubModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ScaleForce Declaration
\*---------------------------------------------------------------------------*/

class ScaleForce
:
    public forceSubModel
{
private:

    dictionary propsDict_;

    word scaleFieldName_;

    const volScalarField& scaleField_;

public:

    //- Runtime type information
    TypeName("ScaleForce");


    // Constructors

        //- Construct from components
        ScaleForce
        (
            const dictionary& dict,
            cfdemCloud& sm,
            forceModel& fm
        );

    // Destructor

        ~ScaleForce();


    // Member Functions

        void partToArray(label, vector&, const vector&, const vector& Ufluid=vector::zero, scalar Cd=scalar(0)) const;

        word myType() const{return typeName; };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
