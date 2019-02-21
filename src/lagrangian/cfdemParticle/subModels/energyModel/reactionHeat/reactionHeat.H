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

\*---------------------------------------------------------------------------*/

#ifndef reactionHeat_H
#define reactionHeat_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "energyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class energyModel Declaration
\*---------------------------------------------------------------------------*/

class reactionHeat
:
    public energyModel

{

protected:

    // Protected data

    dictionary propsDict_;

    bool interpolation_;

    bool verbose_;

    const fvMesh& mesh_;

    scalar maxSource_;

    word reactionHeatName_;

    mutable double **reactionHeat_;

    volScalarField reactionHeatField_;

    void allocateMyArrays() const;

public:

    //- Runtime type information
    TypeName("reactionHeat");

    // Constructors

        //- Construct from components
        reactionHeat
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );

    // Destructor

        virtual ~reactionHeat();

    // Member Functions

        void addEnergyContribution(volScalarField&) const;

        void addEnergyCoefficient(volScalarField&) const {}

        void calcEnergyContribution();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
