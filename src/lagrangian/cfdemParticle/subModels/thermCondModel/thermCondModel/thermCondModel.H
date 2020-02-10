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
    model for thermal conductivity of fluid phase in presence of particles
    special treatment for near-wall behavior:
    a) wallQFactor can be used to impose scaling factor for thermal conductivity at the boundaries,
    e.g. to model a thin boundary layer below grid resolution
    b) explicitly specify a boundary layer thickness
    c) specify given heat loss in W / (m K)
    ATTENTION: if zeroGradient BC are used in thCondField, they cannot be overwritten with numerical data

Class
    thermCondModel

SourceFiles
    thermCondModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermCondModel_H
#define thermCondModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class thermCondModel Declaration
\*---------------------------------------------------------------------------*/

class thermCondModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

        const fvMesh& mesh_;

        IOdictionary transportProperties_;

        dimensionedScalar kf0_;

        volScalarField& thermCondField_;

        // ratio of half-cell-size and near-wall film
        mutable volScalarField wallQFactor_;

        bool hasWallQFactor_;

        // thickness of near-wall boundary layer, can be used instead of Q factor
        mutable volScalarField wallBoundaryLayerThickness_;

        bool hasWallBoundaryLayerThickness_;

        // thickness of near-wall boundary layer, can be used instead of Q factor
        mutable volScalarField wallHeatLoss_;

        bool hasWallHeatLoss_;

public:

    //- Runtime type information
    TypeName("thermCondModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermCondModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        thermCondModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~thermCondModel();


    // Selector

        static autoPtr<thermCondModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Member Functions

        virtual void calcThermCond() {}

        virtual void calcBoundaryCorrections();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
