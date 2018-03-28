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

    Copyright (C) 2018- Mathias Vångö, JKU Linz, Austria

Class
    gradPForceSmooth

Description
    This code is an extention of the gradPForce class by allowing the pressure-
    field to be smoothened prior to the force calculation (without altering the
    original field).

SourceFiles
    gradPForceSmooth.C

\*---------------------------------------------------------------------------*/

#ifndef gradPForceSmooth_H
#define gradPForceSmooth_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class gradPForceSmooth Declaration
\*---------------------------------------------------------------------------*/

class gradPForceSmooth
:
    public forceModel
{
private:
    dictionary propsDict_;

    word pFieldName_;

    const volScalarField& p_;

    word velocityFieldName_;

    const volVectorField& U_;

    bool useRho_;

    bool useU_;          // if false: substitution p=0.5*rho*U^2

    mutable double addedMassCoeff_; //added mass coefficient
    
    autoPtr<smoothingModel> smoothingModel_;
    
    mutable volScalarField pSmooth_;    

public:

    //- Runtime type information
    TypeName("gradPForceSmooth");


    // Constructors

        //- Construct from components
        gradPForceSmooth
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~gradPForceSmooth();


    // Member Functions
        void setForce() const;
        inline const smoothingModel& smoothingM() const
        {
            return smoothingModel_;
        }        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
