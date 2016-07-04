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
    fines influence on particles

SourceFiles
    FinesDynFanning.C
\*---------------------------------------------------------------------------*/

#ifndef FinesDynFanning_H
#define FinesDynFanning_H

#include "forceModel.H"
#include "interpolationCellPoint.H"
#include "FinesDynFanningFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FinesDynFanning Declaration
\*---------------------------------------------------------------------------*/

class FinesDynFanning
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
    
    const volVectorField& UDyn_;
    
    const volScalarField& alphaDyn_;
    
    const volScalarField& alphaP_;
    
    const volScalarField& dHydMix_;
    
    const volScalarField& dSauter_;
    
    const volScalarField& Froude_;
    
    scalar rhoDyn_;
    
    scalar prefactor_;
    
    scalar exponent_;

    mutable scalar scaleDia_;

    mutable scalar scaleDrag_;
    

    

public:

    //- Runtime type information
    TypeName("FinesDynFanning");


    // Constructors

        //- Construct from components
        FinesDynFanning
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~FinesDynFanning();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
