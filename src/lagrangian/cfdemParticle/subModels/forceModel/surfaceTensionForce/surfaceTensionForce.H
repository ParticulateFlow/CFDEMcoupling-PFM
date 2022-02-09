/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code. If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2018- Mathias Vångö, JKU Linz, Austria

Description
    Applies a pre-calculated surface tension force to the particles.

Class
    surfaceTensionForce

SourceFiles
    surfaceTensionForce.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceTensionForce_H
#define surfaceTensionForce_H

#include "forceModel.H"
#include "interpolationCellPoint.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class surfaceTensionForce Declaration
\*---------------------------------------------------------------------------*/

class surfaceTensionForce
:
    public forceModel
{
private:
    dictionary propsDict_;

    word stfFieldName_;
    
    const surfaceScalarField& stf_;

public:

    //- Runtime type information
    TypeName("surfaceTensionForce");


    // Constructors

        //- Construct from components
        surfaceTensionForce
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~surfaceTensionForce();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
