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

#ifndef freeStreaming_H
#define freeStreaming_H

#include "forceModelRec.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class freeStreaming Declaration
\*---------------------------------------------------------------------------*/

class freeStreaming
:
    public forceModelRec
{
private:
    dictionary propsDict_;
    
    bool interpolate_;
    
    word UsFieldName_;
    
    const volVectorField& Us_;
    
    word UsRecFieldName_;
    
    const volVectorField& UsRec_;
    
    word voidfractionRecFieldName_;
    
    const volScalarField& voidfractionRec_;
    
    scalar critVoidfraction_;
    
    scalar particleDensity_;
    
    const vector gravAcc_;

public:

    //- Runtime type information
    TypeName("freeStreaming");


    // Constructors

        //- Construct from components
        freeStreaming
        (
            const dictionary& dict,
            cfdemCloudRec& sm
        );

    // Destructor

        ~freeStreaming();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
