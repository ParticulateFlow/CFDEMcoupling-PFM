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
#include "multiIntervalPath.H"
#include "Random.H"
#include "recModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(multiIntervalPath, 0);

addToRunTimeSelectionTable
(
    recPath,
    multiIntervalPath,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
multiIntervalPath::multiIntervalPath
(
    const dictionary& dict,
    recBase& base
)
:
    recPath(dict, base),
    propsDict_(dict.subDict(typeName + "Props")),
    meanIntervalSteps_(propsDict_.lookupOrDefault<label>("meanIntervalSteps",-1)),
    numIntervals_(0),
    intervalSizes_(propsDict_.lookup("intervalSizes")),
    intervalSizesCumulative_(intervalSizes_),
    Pjump_(0.0),
    intervalWeights_(propsDict_.lookup("intervalWeights")),
    intervalWeightsCumulative_(intervalWeights_)
{
    numIntervals_ = intervalWeights_.size();
    
    if(numIntervals_ == 0 || numIntervals_ != intervalSizes_.size())
    {
        FatalError << ": bad number of separation times and/or interval weights!\n"
                       << abort(FatalError);
    }
    
    if(meanIntervalSteps_<0)
    {
        // if no mean interval length for consecutive steps is specified, use 1/5 from first interval
        meanIntervalSteps_ = (label) (0.2 * intervalSizes_[0]);
    }
    
    // normalize weights
    scalar wsum = 0.0;
    for(int i=0;i<numIntervals_;i++)
    {
        wsum += intervalWeights_[i];
    }
    
    for(int i=0;i<numIntervals_;i++)
    {
        intervalWeights_[i] /= wsum;
    }
    
    for(int i=0;i<numIntervals_;i++)
    {
        scalar sum1 = 0.0;
	scalar sum2 = 0.0;
	for(int j=0;j<=i;j++)
	{
	    sum1 += intervalWeights_[j];
	    sum2 += intervalSizes_[j];
	}
	intervalWeightsCumulative_[i] = sum1;
	intervalSizesCumulative_[i] = sum2;
    }
    
    label numRecFields = base_.recM().numRecFields();
    if (numRecFields != intervalSizesCumulative_[numIntervals_ - 1])
    {
      FatalError << "total number of recurrence fields " << numRecFields << " does not match sum of intervals "
          << intervalSizesCumulative_[numIntervals_ - 1] << abort(FatalError);
    }
    
    // given a jump probability of P, the probability of finding a chain of length N is
    // P(N) = (1 - P)^N * P, and the mean length E(N) = (1 - P) / P
    Pjump_ = 1.0 / (1 + meanIntervalSteps_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

multiIntervalPath::~multiIntervalPath()
{}

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * //

void multiIntervalPath::getRecPath()
{
    label numRecIntervals = 0;
    
    if (Pstream::master())
    {
        computeRecPath();
        numRecIntervals=virtualTimeIndexList_.size();
    }
    
    Pstream::scatter(numRecIntervals);
    
    if (not Pstream::master())
    {
        virtualTimeIndexList_.setSize(numRecIntervals);
    }
    
    Pstream::scatter(virtualTimeIndexList_);
    
    if (verbose_)
    {
        Info << "\nRecurrence path communicated to all processors.\n" << endl;
    }
}


void multiIntervalPath::computeRecPath()
{
    Info << "\nComputing recurrence path\n" << endl;
    
    
    Random ranGen(osRandomInteger());
    
    label virtualTimeIndex=0;
    label recSteps=0;
    label seqStart=0;
    bool prevStepWasJump = true; 

    SymmetricSquareMatrix<scalar>& recurrenceMatrix( base_.recM().recurrenceMatrix() );


    if(base_.recM().totRecSteps() == 1)
    {
        Info<< "\nPrimitive recurrence path with one element.\n" << endl;
        return;
    }
   
    while(recSteps <= base_.recM().totRecSteps() )
    {
	scalar randJump = ranGen.scalar01();
	
	// check if current virtualTimeIndex is close to separation time
	bool intervalBorder = false;
	label sep = 0;
	for(int i = 0;i < numIntervals_; i++)
	{
	    sep += intervalSizes_[i];
	    if (sep - 1 == virtualTimeIndex) intervalBorder=true;
	}
	
	if ((randJump > Pjump_ && !intervalBorder) || prevStepWasJump)
	{
	    virtualTimeIndex++;
	    //recSteps++;
	    prevStepWasJump = false;
	}
	else
	{
	    // before jump, complete former consecutive interval
	    labelPair seqStartEnd(seqStart,virtualTimeIndex);
	    virtualTimeIndexList_.append(seqStartEnd);
	    recSteps += virtualTimeIndex - seqStart + 1;
	    
	    // now jump
	    
	    // identify interval to jump to
	    scalar randInterval = ranGen.scalar01();
	    
	    label interval = numIntervals_-1;
	    for(int i = numIntervals_-2 ;i >= 0; i--)
	    {
	        if (randInterval < intervalWeightsCumulative_[i]) interval=i;
	    }
	    
	    label startLoop = 0;
	    if (interval > 0) startLoop = intervalSizesCumulative_[interval-1];
            label endLoop = intervalSizesCumulative_[interval] - meanIntervalSteps_;
	    
	    scalar nextMinimum(GREAT);
            for (label j = startLoop; j <= endLoop; j++)
	    {
	        if(abs(j - virtualTimeIndex) < meanIntervalSteps_) continue;
		if (recurrenceMatrix[j][virtualTimeIndex] < nextMinimum)
        	{
        		nextMinimum = recurrenceMatrix[j][virtualTimeIndex];
        		seqStart = j+1;
        	}
	    }
	    virtualTimeIndex = seqStart;
	    prevStepWasJump = true;
	}
    }
    Info<< "\nComputing recurrence path done\n" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
