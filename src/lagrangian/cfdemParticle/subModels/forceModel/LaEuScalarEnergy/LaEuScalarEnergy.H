/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
                                Thomas Lichtenegger, thomas.lichtenegger@jku.at
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-2015 DCS Computing GmbH, Linz
                                Copyright 2015-     JKU Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

    two way LaEu Scalar Exchange between DEM and CFD
    convective heat transfer model following
    Li and Mason (2000), A computational investigation of transient heat
    transfer in  pneumatic transport of granular particles, Pow.Tech 112

Class
    LaEuScalarEnergy

SourceFiles
    LaEuScalarEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef LaEuScalarEnergy_H
#define LaEuScalarEnergy_H

#include "forceModel.H"
#include "averagingModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LaEuScalarEnergy Declaration
\*---------------------------------------------------------------------------*/

class LaEuScalarEnergy
:
    public forceModel
{
private:

    dictionary propsDict_;

    word tempFieldName_;

    const volScalarField& tempField_;            // ref to temperature field

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;    // ref to voidfraction field

    scalar maxSource_;                           // max (limited) value of src field

    word velFieldName_;

    const volVectorField& U_;

    word partTempName_;

    mutable double **partTemp_;                          // Lagrangian array

    word partHeatFluxName_;

    mutable double **partHeatFlux_;              // Lagrangian array

    scalar lambda_;                              // fluid thermal conductivity [W/(m*K)]
    
    scalar lambdaPart_;                          // particle thermal conductivity

    scalar Cp_;                                  // specific heat capacity [W*s/(kg*K)]

    void allocateMyArrays() const;

public:

    //- Runtime type information
    TypeName("LaEuScalarEnergy");

    // Constructors

        //- Construct from components
        LaEuScalarEnergy
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~LaEuScalarEnergy();


    // Member Functions
        void setForce() const;

        void partFluidHeatTrans(volScalarField&) const;
	
        scalar thermCondPart() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
