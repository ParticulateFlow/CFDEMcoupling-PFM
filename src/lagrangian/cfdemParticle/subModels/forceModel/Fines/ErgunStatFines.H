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

Class
    ErgunStatFines

SourceFiles
    ErgunStatFines.C

\*---------------------------------------------------------------------------*/

#ifndef ErgunStatFines_H
#define ErgunStatFines_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ErgunStatFines Declaration
\*---------------------------------------------------------------------------*/

class ErgunStatFines
:
    public forceModel
{
private:
    dictionary propsDict_;

    word velFieldName_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    const volScalarField& dSauter_;

    const volScalarField& dSauterMix_;

    const volScalarField& alphaP_;

    const volScalarField& alphaSt_;

    const scalar phi_;

    word UsFieldName_;

    const volVectorField& UsField_;  // the average particle velocity field

    mutable scalar scaleDia_;

    mutable scalar scaleDist_;

    mutable scalar scaleDrag_;

    mutable scalar switchingVoidfraction_; //voidfraction above which dilute formulation will be used

    inline scalar dSauter(label) const;

public:

    //- Runtime type information
    TypeName("ErgunStatFines");


    // Constructors

        //- Construct from components
        ErgunStatFines
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~ErgunStatFines();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
