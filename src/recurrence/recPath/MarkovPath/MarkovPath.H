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

Description
    A recurrence database consisting of N separate intervals is assumed with separation
    times t0 (start time), t1, ... tN (end time) and weights w0, ... wN-1.
\*---------------------------------------------------------------------------*/

#ifndef MarkovPath_H
#define MarkovPath_H

#include "recPath.H"
#include "scalarList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class MarkovPath Declaration
\*---------------------------------------------------------------------------*/

class MarkovPath
:
    public recPath
{
protected:

    // Protected data

    dictionary propsDict_;

    void extendPath();

    void getRecPath();

    void computeRecPath();

    void weightsCumulation();

    void weightsNormalization();
    
    // correlation steps determine which temporal neighborhood is excluded from nearest-neighbor search
    label correlationSteps_;

    label meanIntervalSteps_;
    
    label minIntervalSteps_;

    label numIntervals_;

    label recSteps_;

    label startIndex_;

    labelList intervalSizes_;

    labelList intervalSizesCumulative_;

    scalar Pjump_;

    scalarList intervalWeights_;

    scalarList intervalWeightsCumulative_;

    Random ranGen;


public:

    //- Runtime type information
    TypeName("MarkovPath");

    // Constructors

        //- Construct from components
        MarkovPath
        (
            const dictionary& dict,
            recBase& base
        );


    // Destructor

        virtual ~MarkovPath();


    // Member Functions

    void updateIntervalWeights(scalarList);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
