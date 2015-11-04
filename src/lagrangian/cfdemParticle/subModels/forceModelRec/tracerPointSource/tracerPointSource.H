/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling
    
    Contributing authors:
    Thomas Lichtenegger, Gerhard Holzinger
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


Class
    tracerPointSource

SourceFiles
    tracerPointSource.C

\*---------------------------------------------------------------------------*/

#ifndef tracerPointSource_H
#define tracerPointSource_H

#include "forceModelRec.H"
#include "averagingModel.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class tracerPointSource Declaration
\*---------------------------------------------------------------------------*/

class tracerPointSource
:
    public forceModelRec
{
private:
    dictionary propsDict_;
    
    mutable volScalarField tracerConcentrationField_;       // Eulerian field
    
    mutable double **tracerConcentrationPart_;              // Lagrangian array

    vector sourcePosition_;
    
    scalar sourceStrength_;
    
    scalar beginTime_;
    
    scalar endTime_;
    
    label numReceiverPart_;
    
    mutable std::vector<int> receiverIDs_;
    
    mutable std::vector<scalar> receiverDist_;
    
    mutable std::vector<int>::iterator itID_;
    
    mutable std::vector<scalar>::iterator itDist_;

    void allocateMyArrays() const;

public:

    //- Runtime type information
    TypeName("tracerPointSource");

    // Constructors

        //- Construct from components
        tracerPointSource
        (
            const dictionary& dict,
            cfdemCloudRec& sm
        );

    // Destructor

        ~tracerPointSource();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
