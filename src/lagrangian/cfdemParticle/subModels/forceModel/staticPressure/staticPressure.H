/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    Copyright (C) 2023  Behrad Esgandari, JKU Linz, Austria
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    A gravitational force applied on the particles just in case of
    "Periodic Box simulations"

Class
    staticPressure

SourceFiles
    staticPressure.C

\*---------------------------------------------------------------------------*/

#ifndef staticPressure_H
#define staticPressure_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class staticPressure Declaration
\*---------------------------------------------------------------------------*/

class staticPressure
:
    public forceModel
{
private:

    dictionary propsDict_;

    word rhoGFieldName_;

    scalar rhoGas_;

    vector g_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    scalar solidFraction_;

public:

    //- Runtime type information
    TypeName("staticPressure");


    // Constructors

        //- Construct from components
        staticPressure
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~staticPressure();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
