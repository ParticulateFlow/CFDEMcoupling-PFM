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
    weightSecondaryPhase

SourceFiles
    weightSecondaryPhase.C

\*---------------------------------------------------------------------------*/

#ifndef weightSecondaryPhase_H
#define weightSecondaryPhase_H

#include "otherForceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class weightSecondaryPhase Declaration
\*---------------------------------------------------------------------------*/

class weightSecondaryPhase
:
    public otherForceModel
{

protected:

    // Protected data
        dictionary propsDict_;

        word volfracFieldName_;

        const volScalarField& alpha_;

        dimensionedScalar rho_;

        const dimensionedVector g_;


public:

    //- Runtime type information
    TypeName("weightSecondaryPhase");

    // Constructors

        //- Construct from components
        weightSecondaryPhase
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~weightSecondaryPhase();


    // Member Functions


        tmp<volVectorField> exportForceField();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
