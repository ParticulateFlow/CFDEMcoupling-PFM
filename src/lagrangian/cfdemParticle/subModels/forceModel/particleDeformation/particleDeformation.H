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
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

    depending on a given field and thresholds, particles start deforming
    deformation is characterized by a number in [0,1]

Class
    particleDeformation

SourceFiles
    particleDeformation.C

\*---------------------------------------------------------------------------*/

#ifndef particleDeformation_H
#define particleDeformation_H

#include "forceModel.H"
#include "averagingModel.H"
#include "interpolationCellPoint.H"
#include "autoPtr.H"
#include "cellSet.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class particleDeformation Declaration
\*---------------------------------------------------------------------------*/

class particleDeformation
:
    public forceModel
{
private:

    dictionary propsDict_;

    mutable bool initialExec_;

    word refFieldName_;

    mutable autoPtr<volScalarField> refField_;

    // default deformation in region
    word defaultDeformCellsName_;

    autoPtr<cellSet> defaultDeformCells_;

    bool existDefaultDeformCells_;

    scalar defaultDeformation_;

    labelList partTypes_;

    scalarList lowerBounds_;

    scalarList upperBounds_;

    const word partDeformationsName_;

    label getListIndex(label) const;

    void init() const;

    bool defaultDeformCell(label) const;

public:

    //- Runtime type information
    TypeName("particleDeformation");

    // Constructors

        //- Construct from components
        particleDeformation
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~particleDeformation();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
