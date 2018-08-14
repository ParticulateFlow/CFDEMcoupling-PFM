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
    Model to communicate available reactant per particle

\*---------------------------------------------------------------------------*/

#ifndef reactantPerParticle_H
#define reactantPerParticle_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "chemistryModel.H"
#include "diffusionCoefficients.H"

#include "HashPtrTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class reactantPerParticle Declaration
\*---------------------------------------------------------------------------*/

// for future use:
// + communicate every N steps

class reactantPerParticle
:
    public chemistryModel
{
private:

    dictionary propsDict_;

    const fvMesh& mesh_;

    bool verbose_;

    mutable double **reactantPerParticle_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    mutable volScalarField particlesPerCell_;

    void allocateMyArrays() const;

    label loopCounter_;

    label Nevery_;

public:

    //- Runtime type information
    TypeName("reactantPerParticle");

    // Constructors

        //- Construct from components
        reactantPerParticle
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );

    // Destructor

        virtual ~reactantPerParticle();

    // Member Functions
        void execute();

        void reAllocMyArrays() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
