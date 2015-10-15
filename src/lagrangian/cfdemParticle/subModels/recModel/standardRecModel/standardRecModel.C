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
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "Random.H"
#include "standardRecModel.H"
#include "addToRunTimeSelectionTable.H"
#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardRecModel, 0);

addToRunTimeSelectionTable
(
    recModel,
    standardRecModel,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
standardRecModel::standardRecModel
(
    const dictionary& dict,
    cfdemCloudRec& sm
)
:
    recModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    UFieldName_(propsDict_.lookup("velRecFieldName")),
    voidfractionFieldName_(propsDict_.lookup("voidfractionRecFieldName")),
    UsFieldName_(propsDict_.lookup("granVelRecFieldName")),
    voidfractionRecpl(numRecFields),
    URecpl(numRecFields),
    UsRecpl(numRecFields),
    voidfractionRec_(NULL),
    URec_(NULL),
    UsRec_(NULL)
{
    readFieldSeries();
    
    // make sure each processor has the same sequence of fields
    int root=0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    computeRecMatrix();

    if(verbose_)
      Info << "\nSumming recurrence matrix over all processors.\n" << endl;
    
    for(int i = 0; i < numRecFields; ++i)
      MPI_Allreduce
      (
        &recurrenceMatrixLocal[i][0],
        &recurrenceMatrix[i][0],
        numRecFields,
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD
      );
    
    // root proc computes recurrence path, other processors need to reserve enough space
    if(rank==root)
    {
      computeRecPath();
      numRecIntervals=virtualTimeIndexList.size();
    }
    
    MPI_Bcast(&numRecIntervals,sizeof(numRecIntervals),MPI_BYTE,root,MPI_COMM_WORLD);
    
    if(rank!=root)
      virtualTimeIndexList.setSize(numRecIntervals);
    
    int listSizeBytes = 2*sizeof(sequenceStart)*numRecIntervals;
    
    MPI_Bcast(&virtualTimeIndexList[0], listSizeBytes, MPI_BYTE, root, MPI_COMM_WORLD); 
    if(verbose_)
      Info << "\nRecurrence path communicated to all processors.\n" << endl;
    
    sequenceStart=virtualTimeIndexList[0].first();
    sequenceEnd=virtualTimeIndexList[0].second();
    virtualTimeIndex=sequenceStart;
    
    voidfractionRec_.reset(voidfractionRecpl(virtualTimeIndex));
    URec_.reset(URecpl(virtualTimeIndex));
    UsRec_.reset(UsRecpl(virtualTimeIndex));
    
    writeRecMatrix();    
    writeRecPath();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardRecModel::~standardRecModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void standardRecModel::initRecFields()
{
}

void standardRecModel::updateRecFields()
{
  virtualTimeIndex++;
  if(virtualTimeIndex>sequenceEnd)
  {
    virtualTimeIndexListPos++;
    sequenceStart=virtualTimeIndexList[virtualTimeIndexListPos].first();
    sequenceEnd=virtualTimeIndexList[virtualTimeIndexListPos].second();
    virtualTimeIndex=sequenceStart;
  }
  voidfractionRec_.reset(voidfractionRecpl(virtualTimeIndex));
  URec_.reset(URecpl(virtualTimeIndex));
  UsRec_.reset(UsRecpl(virtualTimeIndex));
}

void standardRecModel::writeRecFields() const
{
 // need to check if this writes to the correct destination 
  
 // voidfractionRecpl[virtualTimeIndex].writeRecFields();
 // URecpl[virtualTimeIndex].writeRecFields();
 // UsRecpl[virtualTimeIndex].writeRecFields();
}

autoPtr<const volScalarField> standardRecModel::voidfraction() const
{
  return voidfractionRec_;
}

autoPtr<const volVectorField> standardRecModel::U() const
{
  return URec_;
}

autoPtr<const volVectorField> standardRecModel::Us() const
{
  return UsRec_;
}

tmp<volScalarField> standardRecModel::tvoidfraction() const
{
  return voidfractionRecpl[virtualTimeIndex];
}

tmp<volVectorField> standardRecModel::tU() const
{
  return URecpl[virtualTimeIndex];
}

tmp<volVectorField> standardRecModel::tUs() const
{
  return UsRecpl[virtualTimeIndex];
}

void standardRecModel::readFieldSeries()
{
     Info << "Reading fields" << endl;
    
    for (instantList::iterator it=timeDirs.begin(); it != timeDirs.end(); ++it)
    {
        // set time
        recTime.setTime(*it, it->value());
        
        // skip constant
        if (recTime.timeName() == "constant")
        {
        	continue;
        }
        
        if (verbose_)
    	{
        	Info << "Reading at t = " << recTime.timeName() << endl;
        }
        
        voidfractionRecpl.set
        (
            timeIndexList(recTime.timeName()),
            new volScalarField
            (
                IOobject
                (
                    voidfractionFieldName_,
                    recTime.timePath(),
                    particleCloud_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                particleCloud_.mesh()
            )
        );
            
        URecpl.set
        (
            timeIndexList(recTime.timeName()),
            new volVectorField
            (
                IOobject
                (
                    UFieldName_,
                    recTime.timePath(),
                    particleCloud_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                particleCloud_.mesh()
            )
        );
        UsRecpl.set
        (
            timeIndexList(recTime.timeName()),
            new volVectorField
            (
                IOobject
                (
                    UsFieldName_,
                    recTime.timePath(),
                    particleCloud_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                particleCloud_.mesh()
            )
        );
        
    }
    Info << "Reading fields done" << endl;
}

void standardRecModel::computeRecMatrix()
{
    Info<< "\nComputing recurrence matrix\n" << endl;
    //scalar maxElemVal(0.0);
    
    // compute recurrence matrix elements
    forAll(timeIndexList, ti)
    {
    	forAll(timeIndexList, tj)
    	{
	        if(verbose_)
		  Info<<"\n Doing calculation for element " << ti << " " << tj << "\n" << endl;
    		// main diagonal
    		if (ti == tj)
    		{
    			recurrenceMatrixLocal[ti][tj] = 0;
    			continue;
    		}
    		
    		// skip one half of the matrix
    		if (ti > tj)
    		{
    			continue;
    		}
    		
    		// compute elements
    		recurrenceMatrixLocal[ti][tj]
    			= sumSqr(voidfractionRecpl[ti].internalField() - voidfractionRecpl[tj].internalField());
    		
    		recurrenceMatrixLocal[tj][ti] = recurrenceMatrixLocal[ti][tj];
    		
    	//	if (maxElemVal < recurrenceMatrixLocal[ti][tj])
    	//	{
    	//		maxElemVal = recurrenceMatrixLocal[ti][tj];
    	//	}
    	}
    }
    
    // normalize matrix elements
    //forAll(timeIndexList, ti)
    //{
    //	forAll(timeIndexList, tj)
    //	{
    //		recurrenceMatrixLocal[ti][tj] /= maxElemVal;
    //	}
    //} 
    Info<< "\nComputing recurrence matrix done\n" << endl;
}

void standardRecModel::computeRecPath()
{
    Info<< "\nComputing recurrence path\n" << endl;
    Random ranGen(osRandomInteger());
    
    label virtualTimeIndex=0;
    label recSteps=0;
    label seqStart=0;
    label seqLength=ranGen.integer(lowerSeqLim, upperSeqLim);
    virtualTimeIndex=seqEnd(seqStart,seqLength);
    labelPair seqStartEnd(seqStart,virtualTimeIndex);
    virtualTimeIndexList.append(seqStartEnd);
    recSteps+=seqLength;
   
    while(recSteps<=totRecSteps)
    {
        label startLoop = 0;
        label endLoop = 0;

        // search the other half of the recurrence matrix for 
        // the new starting point of the next sequence
        if (virtualTimeIndex < recurrenceMatrix.n()/2)
        {
        	startLoop = recurrenceMatrix.n()/2;
        	endLoop = recurrenceMatrix.n()-1;
		endLoop--; // start of next sequence one snapshot AFTER minimum position
        }
        else
        {
        	startLoop = 0;
        	endLoop = recurrenceMatrix.n()/2-1;
        }
        

        scalar nextMinimum(GREAT);
        for (label j = startLoop; j < endLoop; j++)
        {
        	if (recurrenceMatrix[j][virtualTimeIndex] < nextMinimum)
        	{
        		nextMinimum = recurrenceMatrix[j][virtualTimeIndex];
        		seqStart = j+1;
        		continue;
        	}
        }
        
        seqLength = ranGen.integer(lowerSeqLim, upperSeqLim);  
        virtualTimeIndex=seqEnd(seqStart,seqLength);
        labelPair seqStartEnd(seqStart,virtualTimeIndex);
        virtualTimeIndexList.append(seqStartEnd);
        recSteps+=seqLength;
    }
    Info<< "\nComputing recurrence path done\n" << endl;
}

label standardRecModel::seqEnd(label seqStart, label & seqLength)
{
  if(seqStart+seqLength>numRecFields-1)
      seqLength=numRecFields-1-seqStart;
  return seqStart+seqLength;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
