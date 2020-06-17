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

#ifndef isotropicFluctuations_H
#define isotropicFluctuations_H

#include "forceModel.H"
#include "interpolationCellPoint.H"
#include "Random.H"
#include "autoPtr.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class isotropicFluctuations Declaration
\*---------------------------------------------------------------------------*/

class isotropicFluctuations
:
    public forceModel
{
protected:
    dictionary propsDict_;

    bool interpolate_;

    bool measureDiff_;

    mutable OFstream recErrorFile_;

    // ignore fluctuations in region
    word ignoreCellsName_;

    autoPtr<cellSet> ignoreCells_;

    bool existIgnoreCells_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    word voidfractionRecFieldName_;

    const volScalarField& voidfractionRec_;

    scalar critVoidfraction_;

    scalar D0_;

    word D0FieldName_;

    volScalarField D0Field_;

    scalar maxDisplacement_;

    const scalar dtDEM_;

    virtual vector unitFlucDir() const;

    scalar fluctuationMag(const scalar, const scalar) const;

    mutable Random ranGen_;

    bool ignoreCell(label) const;

public:

    //- Runtime type information
    TypeName("isotropicFluctuations");


    // Constructors

        //- Construct from components
        isotropicFluctuations
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word type = "isotropicFluctuations"
        );

    // Destructor

        ~isotropicFluctuations();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
