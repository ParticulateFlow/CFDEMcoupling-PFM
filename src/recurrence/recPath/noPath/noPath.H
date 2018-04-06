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

#ifndef noPath_H
#define noPath_H

#include "recPath.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noPath Declaration
\*---------------------------------------------------------------------------*/

class noPath
:
    public recPath
{
protected:

    // Protected data

    dictionary propsDict_;

    void getRecPath();

    void computeRecPath();

    void computeJumpVector();

    label seqEnd(label, label&);



public:

    //- Runtime type information
    TypeName("noPath");

    // Constructors

        //- Construct from components
        noPath
        (
            const dictionary& dict,
            recBase& base
        );


    // Destructor

        virtual ~noPath();


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
