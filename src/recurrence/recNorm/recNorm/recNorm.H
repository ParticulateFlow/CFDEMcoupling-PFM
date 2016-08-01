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

#ifndef recNorm_H
#define recNorm_H

#include "recBase.H"
#include "fvCFD.H"
#include "HashTable.H"
#include "labelPair.H"
#include "OFstream.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class recNorm Declaration
\*---------------------------------------------------------------------------*/

class recNorm
{
protected:

    // Protected data
    recBase& base_;
    
    const dictionary& recProperties_;
    
    IOdictionary controlDict_;
    
    bool verbose_;
    

    
   

public:

    //- Runtime type information
    TypeName("recNorm");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            recNorm,
            dictionary,
            (
                const dictionary& dict,
                recBase& base
            ),
            (dict,base)
        );


    // Constructors

        //- Construct from components
        recNorm
        (
            const dictionary& dict,
            recBase& base
        );


    // Destructor

        virtual ~recNorm();


    // Selector

        static autoPtr<recNorm> New
        (
            const dictionary& dict,
            recBase& base
        );


    // Member Functions
    
  
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
