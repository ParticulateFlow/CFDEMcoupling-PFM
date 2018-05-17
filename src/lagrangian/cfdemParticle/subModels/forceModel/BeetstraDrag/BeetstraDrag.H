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
    drag law for monodisperse systems according to
    Beetstra et al. AIChE J 53.2 (2007)

SourceFiles
    BeetstraDrag.C
\*---------------------------------------------------------------------------*/

#ifndef BeetstraDrag_H
#define BeetstraDrag_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BeetstraDrag Declaration
\*---------------------------------------------------------------------------*/

class BeetstraDrag
:
    public forceModel
{
private:
    dictionary propsDict_;

    word velFieldName_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    word UsFieldName_;

    const volVectorField& UsField_;

    mutable scalar scaleDia_;

    mutable scalar scaleDrag_;

    mutable bool polydisperse_;

    word dSauterFieldName_;

    const volScalarField& dSauterField_;

public:

    //- Runtime type information
    TypeName("BeetstraDrag");


    // Constructors

        //- Construct from components
        BeetstraDrag
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~BeetstraDrag();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
