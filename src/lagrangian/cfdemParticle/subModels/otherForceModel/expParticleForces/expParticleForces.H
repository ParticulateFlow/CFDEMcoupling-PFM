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
    expParticleForces

SourceFiles
    expParticleForces.C

\*---------------------------------------------------------------------------*/

#ifndef expParticleForces_H
#define expParticleForces_H

#include "otherForceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class expParticleForces Declaration
\*---------------------------------------------------------------------------*/

class expParticleForces
:
    public otherForceModel
{

protected:

    // Protected data
    //    dictionary propsDict_;

public:

    //- Runtime type information
    TypeName("expParticleForces");

    // Constructors

        //- Construct from components
        expParticleForces
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~expParticleForces();



    // Member Functions


        tmp<volVectorField> exportForceField();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
