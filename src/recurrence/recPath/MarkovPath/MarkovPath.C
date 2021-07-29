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

#include "error.H"
#include "MarkovPath.H"
#include "Random.H"
#include "recModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(MarkovPath, 0);

addToRunTimeSelectionTable
(
    recPath,
    MarkovPath,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
MarkovPath::MarkovPath
(
    const dictionary& dict,
    recBase& base
)
:
    recPath(dict, base),
    propsDict_(dict.subDict(typeName + "Props")),
    searchMinimum_(propsDict_.lookupOrDefault<bool>("searchMinimum",true)),
    correlationSteps_(readLabel(propsDict_.lookup("correlationSteps"))),
    meanIntervalSteps_(propsDict_.lookupOrDefault<label>("meanIntervalSteps",-1)),
    minIntervalSteps_(propsDict_.lookupOrDefault<label>("minIntervalSteps",0)),
    numIntervals_(base.recM().numIntervals()),
    recSteps_(0),
    startIndex_(propsDict_.lookupOrDefault<label>("startIndex",0)),
    intervalSizes_(numIntervals_),
    intervalSizesCumulative_(numIntervals_),
    Pjump_(0.0),
    intervalWeights_(propsDict_.lookupOrDefault<scalarList>("intervalWeights",scalarList(numIntervals_,1.0))),
    intervalWeightsCumulative_(intervalWeights_),
    ranGen(clock::getTime()+pid())
{
    for(int i=0;i<numIntervals_;i++)
    {
        intervalSizes_[i] = base.recM().numRecFields(i);
    }

    if(meanIntervalSteps_<0)
    {
        // if no mean interval length for consecutive steps is specified, use 1/5 from first interval
        meanIntervalSteps_ = static_cast<label>(0.2 * intervalSizes_[0]);
    }

    // normalize weights and create cumulative distribution
    weightsNormalization();
    weightsCumulation();

    for(int i=0;i<numIntervals_;i++)
    {
        scalar sum1 = 0.0;
        for(int j=0;j<=i;j++)
        {
            sum1 += intervalSizes_[j];
        }
        intervalSizesCumulative_[i] = sum1;
    }

    // check if meanIntervalSteps and correlationSteps are reasonable
    label critLength = meanIntervalSteps_ + 2 * correlationSteps_;
    for(int i=0;i<numIntervals_;i++)
    {
        if (critLength >= intervalSizes_[i])
        {
            FatalError <<"too big mean interval size and correlation time for database " << i << "\n" << abort(FatalError);
        }
    }

    // given a jump probability of P, the probability of finding a chain of length N is
    // P(N) = (1 - P)^N * P, and the mean length E(N) = (1 - P) / P
    Pjump_ = 1.0 / (1 + meanIntervalSteps_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

MarkovPath::~MarkovPath()
{}

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * //

void MarkovPath::getRecPath()
{
    recPath::getRecPath();
    Pstream::scatter(recSteps_);
}


void MarkovPath::computeRecPath()
{
    if(base_.recM().totRecSteps() == 1)
    {
        Info << "\nPrimitive recurrence path with one element.\n" << endl;
        return;
    }

    // extend path if
    // a) it does already exist --> add a single interval
    // b) it does not exist yet --> until number of initial steps is reached
    if(recSteps_ > base_.recM().totRecSteps() )
    {
        extendPath();
        base_.recM().writeRecPathLastInterval();
        Info << "\nExtending recurrence path done\n" << endl;
        return;
    }

    while(recSteps_ <= base_.recM().totRecSteps() )
    {
        extendPath();
    }
    Info << "\nComputing recurrence path done\n" << endl;
}

void MarkovPath::extendPath()
{

    SymmetricSquareMatrix<scalar>& recurrenceMatrix( base_.recM().recurrenceMatrix() );

    label seqStart=0;
    label seqEnd=0;
    label virtualTimeIndex=0;

    // if previous intervals exist, perform a jump, otherwise start with step 0
    if(virtualTimeIndexList_.size() > 0 )
    {
        virtualTimeIndex = virtualTimeIndexList_.last().second();

        // jump to similar state in same or other database
        scalar randInterval = ranGen.scalar01();

        label interval = numIntervals_-1;
        for(int i = numIntervals_-2; i >= 0; i--)
        {
           if (randInterval < intervalWeightsCumulative_[i]) interval=i;
        }

        label startLoop = 0;
        if (interval > 0) startLoop = intervalSizesCumulative_[interval-1];
        label endLoop = intervalSizesCumulative_[interval] - meanIntervalSteps_;

        if (searchMinimum_)
        {
            scalar nextMinimum(GREAT);
            for (label j = startLoop; j <= endLoop; j++)
            {
                if(abs(j - virtualTimeIndex) < correlationSteps_) continue;
                if (recurrenceMatrix[j][virtualTimeIndex] < nextMinimum)
                {
                    nextMinimum = recurrenceMatrix[j][virtualTimeIndex];
                    seqStart = j+1;
                }
            }
        }
        else
        {
            scalar randIntStart = ranGen.scalar01();
            seqStart = static_cast<label>(startLoop + randIntStart*(endLoop - startLoop));
        }

        virtualTimeIndex = seqStart;
    }
    else
    {
        seqStart=startIndex_;
    }

    // take a series of consecutive steps
    bool takeAnotherStep = true;
    while(takeAnotherStep)
    {
        virtualTimeIndex++;

        // take another step according to jump probability? only if minIntervalSteps taken
        scalar randJump = ranGen.scalar01();
        if (randJump < Pjump_ && virtualTimeIndex - seqStart >= minIntervalSteps_) takeAnotherStep=false;

        // interval border? must jump
        for(int i = 0;i < numIntervals_; i++)
        {
            if (intervalSizesCumulative_[i] - 1 == virtualTimeIndex) takeAnotherStep=false;
        }
    }

    seqEnd = virtualTimeIndex;

    // add interval to recurrence path
    labelPair seqStartEnd(seqStart,seqEnd);
    virtualTimeIndexList_.append(seqStartEnd);
    recSteps_ += seqEnd - seqStart + 1;
}


void MarkovPath::weightsCumulation()
{
    for(int i=0;i<numIntervals_;i++)
    {
        scalar sum1 = 0.0;
        for(int j=0;j<=i;j++)
        {
            sum1 += intervalWeights_[j];
        }
        intervalWeightsCumulative_[i] = sum1;
    }
}


void MarkovPath::weightsNormalization()
{
    scalar wsum = 0.0;
    for(int i=0;i<numIntervals_;i++)
    {
        wsum += intervalWeights_[i];
    }

    for(int i=0;i<numIntervals_;i++)
    {
        intervalWeights_[i] /= wsum;
    }
}


void MarkovPath::updateIntervalWeights(scalarList newWeights)
{
    // check if number of weights is correct
    if(newWeights.size() > numIntervals_)
    {
        FatalError <<"more weights than databases specified\n" << abort(FatalError);
    }

    for(int i=0; i<numIntervals_; i++)
    {
        if (i < newWeights.size()) intervalWeights_[i] = newWeights[i];
        else intervalWeights_[i] = 0.0;
    }

    weightsNormalization();
    weightsCumulation();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
