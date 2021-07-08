/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling

    Contributing authors:
    Thomas Lichtenegger
    Copyright (C) 2015- Johannes Kepler University, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling academic.

    CFDEMcoupling academic is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling academic is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling academic.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef turbulentDispersion_H
#define turbulentDispersion_H

#include "forceModel.H"
#include "interpolationCellPoint.H"
#include "Random.H"
#include "autoPtr.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class turbulentDispersion Declaration
\*---------------------------------------------------------------------------*/

class turbulentDispersion
:
    public forceModel
{
protected:
    dictionary propsDict_;

    bool interpolate_;

    const fvMesh& mesh_;

    // ignore fluctuations in region
    word ignoreCellsName_;

    autoPtr<cellSet> ignoreCells_;

    bool existIgnoreCells_;

    bool usePreCalcNut_;

    word nutName_;

    bool usePreCalcK_;

    word kName_;

    bool usePreCalcDispField_;

    word dispFieldName_;

    mutable volVectorField* dispField_;

    mutable volScalarField nut_;

    mutable volScalarField wallIndicatorField_;

    mutable volScalarField delta_;

    scalar Ck_;

    scalar dt_;

    scalar turbulentSchmidtNumber_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    scalar critVoidfraction_;

    virtual vector unitFlucDir() const;

    mutable Random ranGen_;

    bool ignoreCell(label) const;

public:

    //- Runtime type information
    TypeName("turbulentDispersion");


    // Constructors

        //- Construct from components
        turbulentDispersion
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word type = "turbulentDispersion"
        );

    // Destructor

        ~turbulentDispersion();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
