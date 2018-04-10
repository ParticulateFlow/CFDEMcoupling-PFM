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

#ifndef directedDiffusiveRelaxation_H
#define directedDiffusiveRelaxation_H

#include "forceModel.H"
#include "interpolationCellPoint.H"
#include "autoPtr.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class directedDiffusiveRelaxation Declaration
\*---------------------------------------------------------------------------*/

class directedDiffusiveRelaxation
:
    public forceModel
{
private:
    dictionary propsDict_;

    bool interpolate_;

    bool measureDiff_;

    mutable OFstream recErrorFile_;

    // ignore fluctuations in region
    word ignoreCellsName_;

    autoPtr<cellSet> ignoreCells_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    word voidfractionRecFieldName_;

    const volScalarField& voidfractionRec_;

    scalar critVoidfraction_;

    scalar D0_;

    scalar maxDisplacement_;

    mutable volScalarField correctedField_;

    mutable volVectorField relaxStream_;

    mutable volScalarField DField_;

    const scalar dtCFD_;

    const scalar dtDEM_;

    scalar timeFac_;

    scalar relaxForT_;

    void relax(scalar) const;

    bool ignoreCell(label) const;
    

public:

    //- Runtime type information
    TypeName("directedDiffusiveRelaxation");


    // Constructors

        //- Construct from components
        directedDiffusiveRelaxation
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~directedDiffusiveRelaxation();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
