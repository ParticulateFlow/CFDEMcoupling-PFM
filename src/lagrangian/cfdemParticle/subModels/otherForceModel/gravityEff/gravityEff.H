/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    Copyright (C) 2023  Behrad Esgandari, JKU Linz, Austria
-------------------------------------------------------------------------------
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

Description
    Gravitaional source term that is necessary to derive the flow just in case
    of "Periodic Box simulations"

Class
    gravityEff

SourceFiles
    gravityEff.C

\*---------------------------------------------------------------------------*/

#ifndef gravityEff_H
#define gravityEff_H

#include "otherForceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class gravityEff Declaration
\*---------------------------------------------------------------------------*/

class gravityEff
:
    public otherForceModel
{

protected:

    // Protected data
        dictionary propsDict_;

        word voidfractionFieldName_;

        const volScalarField& voidfraction_;

        word rhoGFieldName_;

        const volScalarField& rhoG_;

        word gravityFieldName_;

        const uniformDimensionedVectorField& g_;// ref to gravity

        dimensionedScalar rhoPart_;

        const volScalarField& unity_;
public:

    //- Runtime type information
    TypeName("gravityEff");

    // Constructors

        //- Construct from components
        gravityEff
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor
        virtual ~gravityEff();



    // Member Functions
        tmp<volVectorField> exportForceField();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
