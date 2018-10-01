/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
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

    Koch, Hill drag law
    based on Koch Hill 2001,"Inertial effects in suspensions and porous-media
    flows", Annual Review of fluid mechanics.
    including interpolation of the velocity to the exact position
    including drag coefficient for implicit drag for DEM

Class
    KochHillRWDrag

SourceFiles
    KochHillRWDrag.C

\*---------------------------------------------------------------------------*/

#ifndef KochHillRWDrag_H
#define KochHillRWDrag_H

#include "forceModel.H"
#include "interpolationCellPoint.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class KochHillRWDrag Declaration
\*---------------------------------------------------------------------------*/

class KochHillRWDrag
:
    public forceModel
{
private:
    dictionary propsDict_;

    const bool verbose_;

    word velFieldName_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    word UsFieldName_;

    const volVectorField& UsField_;  // the average particle velocity field

    const bool interpolation_; // use interpolated field values

    mutable scalar scale_;

    scalar cl_;

    scalar rhoP_;

    const bool randomTauE_;

    mutable double **partTime_;                // Lagrangian array

    mutable double **partUfluct_;              // Lagrangian array

    mutable Random RanGen_;

public:

    //- Runtime type information
    TypeName("KochHillRWDrag");


    // Constructors

        //- Construct from components
        KochHillRWDrag
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~KochHillRWDrag();


    // Member Functions
        void setForce() const;

        void reAllocArrays() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
