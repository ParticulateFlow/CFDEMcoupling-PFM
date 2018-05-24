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

#ifndef maxNorm_H
#define maxNorm_H

#include "diffNorm.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class maxNorm Declaration
\*---------------------------------------------------------------------------*/

class maxNorm
:
    public diffNorm
{
protected:

    scalar normVSF(label, label);

    scalar normVVF(label, label);

    scalar normSSF(label, label);

public:

    //- Runtime type information
    TypeName("maxNorm");

    // Constructors

        //- Construct from components
        maxNorm
        (
            const dictionary& dict,
            recBase& base
        );


    // Destructor

        virtual ~maxNorm();


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
