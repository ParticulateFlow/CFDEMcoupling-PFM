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
    gradPStatFines

SourceFiles
    gradPStatFines.C

\*---------------------------------------------------------------------------*/

#ifndef gradPStatFines_H
#define gradPStatFines_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gradPStatFines Declaration
\*---------------------------------------------------------------------------*/

class gradPStatFines
:
    public forceModel
{
private:
    dictionary propsDict_;

    word pFieldName_;

    const volScalarField& p_;

    word velocityFieldName_;

    const volVectorField& U_;
    
    const volScalarField& alphaP_;
    
    const volScalarField& alphaSt_;

    bool useRho_;

    bool useU_;          // if false: substitution p=0.5*rho*U^2

    mutable double addedMassCoeff_; //added mass coefficient

public:

    //- Runtime type information
    TypeName("gradPStatFines");


    // Constructors

        //- Construct from components
        gradPStatFines
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~gradPStatFines();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
