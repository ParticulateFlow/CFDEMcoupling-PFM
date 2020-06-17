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

#ifndef potentialRelaxation_H
#define potentialRelaxation_H

#include "forceModel.H"
#include "interpolationCellPoint.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class potentialRelaxation Declaration
\*---------------------------------------------------------------------------*/

class potentialRelaxation
:
    public forceModel
{
private:
    dictionary propsDict_;
    
    bool interpolate_;
    
    bool measureDiff_;
    
    mutable OFstream recErrorFile_;
    
    word voidfractionFieldName_;
    
    const volScalarField& voidfraction_;
    
    word voidfractionRecFieldName_;
    
    const volScalarField& voidfractionRec_;
    
    scalar critVoidfraction_;
    
    scalar D0_;
    
    scalar D1_;
    
    mutable volScalarField correctedField_;
    
    const scalar dt_;
    
    // ignore particles in cells below plane given with ref point and normal vector
    // normal vector points towards region where fluctuations are permitted
    
    bool ignoreReg_;
    
    vector ignoreDirection_;
    
    vector ignorePoint_;
    
    mutable double **vfluc_;              // Lagrangian array

    void allocateMyArrays() const;
    
    void relax(scalar, scalar) const;

public:

    //- Runtime type information
    TypeName("potentialRelaxation");


    // Constructors

        //- Construct from components
        potentialRelaxation
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~potentialRelaxation();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
