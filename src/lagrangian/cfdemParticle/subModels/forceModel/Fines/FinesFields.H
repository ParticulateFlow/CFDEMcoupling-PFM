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
    class managing fines fields

SourceFiles
    finesFields.C
\*---------------------------------------------------------------------------*/

#ifndef FinesFields_H
#define FinesFields_H

#include "fvCFD.H"
#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FinesFields Declaration
\*---------------------------------------------------------------------------*/

class FinesFields
{
private:
    cfdemCloud& particleCloud_;
  
    dictionary propsDict_;

    word velFieldName_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    word UsFieldName_;   

    const volVectorField& UsField_;
    
    const volScalarField& dSauter_;
    
    volScalarField dSauterMix_;
    
    volScalarField alphaP_;

    volScalarField alphaSt_;
    
    volScalarField alphaDyn_;
    
    volVectorField uDyn_;
    
    volScalarField Sds_;
    
    scalar estimatedVelFrac_;
    
    dimensionedScalar dFine_;
    
    scalar rhoDyn_;
    
    scalar depRate_;
    
    scalar nCrit_;
    
    scalar alphaMax_;
        
    void calcSource();
    
    void updateDSauter();
    
    void updateUDyn();

    void integrateFields();

public:

    //- Runtime type information
    TypeName("FinesFields");


    // Constructors

        //- Construct from components
        FinesFields
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~FinesFields();


    // Member Functions

    
    void update();
    
    scalar alphaDyn(label) const;
    
    scalar rhoDyn() const;
    
    scalar uDyn(label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

