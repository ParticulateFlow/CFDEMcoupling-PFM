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
    transfer fluid properties to LIGGGHTS

SourceFiles
    transferFluidProperties.C
\*---------------------------------------------------------------------------*/

#ifndef transferFluidProperties_H
#define transferFluidProperties_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class transferFluidProperties Declaration
\*---------------------------------------------------------------------------*/

class transferFluidProperties
:
    public forceModel
{
private:

    dictionary propsDict_;

public:

    //- Runtime type information
    TypeName("transferFluidProperties");

    // Constructors

        //- Construct from components
        transferFluidProperties
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~transferFluidProperties();


    // Member Functions
        void setForce() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
