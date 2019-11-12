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

#ifndef readNorm_H
#define readNorm_H

#include "recNorm.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class readNorm Declaration
\*---------------------------------------------------------------------------*/

class readNorm
:
    public recNorm
{
protected:

    // Protected data

    dictionary propsDict_;

    const word recMatName_;

    void computeRecMatrix();


public:

    //- Runtime type information
    TypeName("readNorm");

    // Constructors

        //- Construct from components
        readNorm
        (
            const dictionary& dict,
            recBase& base
        );


    // Destructor

        virtual ~readNorm();


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
