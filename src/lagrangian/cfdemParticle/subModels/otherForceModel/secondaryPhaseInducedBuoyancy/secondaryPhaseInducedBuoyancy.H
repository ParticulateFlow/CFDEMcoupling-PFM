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
    secondaryPhaseInducedBuoyancy

SourceFiles
    secondaryPhaseInducedBuoyancy.C

\*---------------------------------------------------------------------------*/

#ifndef secondaryPhaseInducedBuoyancy_H
#define secondaryPhaseInducedBuoyancy_H

#include "otherForceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class secondaryPhaseInducedBuoyancy Declaration
\*---------------------------------------------------------------------------*/

class secondaryPhaseInducedBuoyancy
:
    public otherForceModel
{

protected:

    // Protected data
        dictionary propsDict_;

        word voidfractionFieldName_;

        const volScalarField& voidfraction_;

        // density dispersed phase - density continuous phase
        dimensionedScalar diffRho_;

        const dimensionedVector g_;


public:

    //- Runtime type information
    TypeName("secondaryPhaseInducedBuoyancy");

    // Constructors

        //- Construct from components
        secondaryPhaseInducedBuoyancy
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~secondaryPhaseInducedBuoyancy();


    // Member Functions


        tmp<volVectorField> exportForceField();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
