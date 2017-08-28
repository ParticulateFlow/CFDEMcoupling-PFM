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
    fines influence on particles

SourceFiles
    Fines.C
\*---------------------------------------------------------------------------*/

#ifndef Fines_H
#define Fines_H

#include "forceModel.H"
#include "interpolationCellPoint.H"
#include "FinesFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Fines Declaration
\*---------------------------------------------------------------------------*/

class Fines
:
    public forceModel
{
private:
    mutable FinesFields finesFields_;

public:

    //- Runtime type information
    TypeName("Fines");


    // Constructors

        //- Construct from components
        Fines
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~Fines();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
