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
                        M. Efe Kinaci, JKU Linz, Austria

    Description
    Model to communicate Sherwood number necessary coefficients

\*---------------------------------------------------------------------------*/

#ifndef massTransferCoeff_H
#define massTransferCoeff_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "chemistryModel.H"
#include "species.H"

#include "HashPtrTable.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class massTransferCoefficient Declaration
\*---------------------------------------------------------------------------*/

class massTransferCoeff
:
    public chemistryModel
{
private:

    dictionary propsDict_;

    bool verbose_;

    bool interpolation_;

    const fvMesh& mesh_;

    word velFieldName_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    word densityFieldName_;

    const volScalarField& rho_;

    word partNuName_;

    mutable double **partNu_;

    word partReynolds_;

    mutable double **partRe_;

    void allocateMyArrays() const;

    mutable scalar scaleDia_;

public:

    //- Runtime type information
    TypeName("massTransferCoeff");

    // Constructors

        //- Construct from components
        massTransferCoeff
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );

    // Destructor

        virtual ~massTransferCoeff();

    // Member Functions
        void execute();

        void reAllocMyArrays() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
