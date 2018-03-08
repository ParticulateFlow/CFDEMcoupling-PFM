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
    recBase& base
)
:
    recModel(dict,base),
    propsDict_(dict.subDict(typeName + "Props")),
    //dataBaseName_("dataBase"),
    dataBaseName_(propsDict_.lookupOrDefault<word>("dataBase", word("dataBase"))),
    //recTime("dataBase", "", "../system", "../constant", false),
    recTime(fileName(dataBaseName_), "", "../system", "../constant", false),
    timeDirs(recTime.times()),
    skipZero_(propsDict_.lookupOrDefault<Switch>("skipZero", Switch(false))),
    numRecFields_(skipZero_ ? label(timeDirs.size())-1 : label(timeDirs.size())),
    recurrenceMatrix_(numRecFields_,scalar(-1.0)),
    timeIndexList_(numRecFields_-1),
    timeValueList_(numRecFields_-1),
    contTimeIndex(0),
    lowerSeqLim_(max(1, label(numRecFields_/20))),
    upperSeqLim_(max(1, label(numRecFields_/5))),
    
    volScalarFieldList_(volScalarFieldNames_.size()),
    volVectorFieldList_(volVectorFieldNames_.size()),
    surfaceScalarFieldList_(surfaceScalarFieldNames_.size())
{   
    if (verbose_)
    {
        // be informative on properties of the "recTime" Time-object
        Info << "recTime.rootPath() " << recTime.rootPath() << endl;
        Info << "recTime.caseName() " << recTime.caseName() << endl;
        Info << "recTime.path() " << recTime.path() << endl;
        Info << "recTime.timePath() " << recTime.timePath() << endl;
        Info << "recTime.timeName() " << recTime.timeName() << endl;
        Info << "timeDirs " << timeDirs << endl;
        Info << "ignoring 0 directory: " << skipZero_ << endl;
    }
    readTimeSeries();

    recTimeStep_ = checkTimeStep();
    totRecSteps_ = 1 + static_cast<label>( (endTime_-startTime_) / recTimeStep_ );

    for(int i=0; i<volScalarFieldNames_.size(); i++)
    {
        volScalarFieldList_[i].setSize(numRecFields_);
    }

    for(int i=0; i<volVectorFieldNames_.size(); i++)
    {
        volVectorFieldList_[i].setSize(numRecFields_);
    }

    for(int i=0; i<surfaceScalarFieldNames_.size(); i++)
    {
        surfaceScalarFieldList_[i].setSize(numRecFields_);
    }

    //   setRecFields();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardRecModel::~standardRecModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
scalar standardRecModel::checkTimeStep()
{
   // check time step of provided data
    scalar dtCur(0.0);
    scalar dtOld(0.0);
    
    if (verbose_)
    {
    	Info << "timeValueList : " << timeValueList_ << endl;
    }
    
    if (timeIndexList_.size() == 1)
    {
        return 1.e10;
    }
    
    forAll(timeValueList_, i)
    {
    	// skip zero
    	if (skipZero_ and timeDirs[i].value() == 0)
    	{
    	    if (verbose_)
    	    {
    	        Info << " ... skipping 0 in checkTimeStep()" << endl;
    	    }
    	    continue;
    	}
    	
    	// compute time step
    	if (timeDirs[i].value() == timeDirs.last().value())
    	{
    		if (verbose_)
    		{
    			Info << ".. leaving loop at " << timeDirs[i] << endl;
    		}
    		// leave loop
    		break;
    	}
    	
    	if (verbose_)
    	{
    		Info << "timeDirs.fcIndex(i)].value(),  timeDirs[i].value() : " 
    			<< timeDirs[timeDirs.fcIndex(i)].value() << "   " << timeDirs[i].value()
    			<< endl;
    	}
    	
    	// the documentation is in the code ;-)
    	//	fcIndex() - return forward circular index, i.e. the next index
        dtCur = timeDirs[timeDirs.fcIndex(i)].value() - timeDirs[i].value();
        
        if (dtOld < SMALL)
        {
        	dtOld = dtCur;
        }
        
        if (abs(dtOld - dtCur) > SMALL)
        {
        	Info << "dtCur, dtOld = " << dtCur << "   " << dtOld << endl;
        	FatalError << "    in setting up data" << nl
				<< "    non-constant time-step of provided simulation data" 
				<< abort(FatalError);
        }
    }
        
    // set deltaT
    recTime.setDeltaT(dtCur, false);
	
	if (verbose_)
    {
		Info << "Setting deltaRecT to " << dtCur << endl;
		Info << "Actual recTime.deltaT = " << recTime.deltaTValue() << endl;
		Info << "Actual runTime.deltaT = " << timeStep_ << endl;
    }
    return dtCur;
}

void standardRecModel::readFieldSeries()
{
    Info << "Reading fields\n" << endl;
    
    label size = timeDirs.size();
    label counter = 0;
    label percentage = 0;
    
    for (instantList::iterator it=timeDirs.begin(); it != timeDirs.end(); ++it)
    {
        if(counter >= 0.1 * percentage * size)
	    {
	        Info << "\t" << 10 * percentage << " \% done" << endl;
	        percentage++;
	    }
	    counter++;
      
        // set time
        recTime.setTime(*it, it->value());
        
        // skip zero
        if (skipZero_ and recTime.timeName() == "0")
        {
            if (verbose_)
    	    {
    	        Info << " ... skipping 0 in readFieldSeries()" << endl;
    	    }
    	    
    	    continue;
        }
        
        // skip constant
        if (recTime.timeName() == "constant")
        {
        	continue;
        }
        
        if (verbose_)
    	{
        	Info << "Reading at t = " << recTime.timeName() << endl;
        }
        
        for(int i=0; i<volScalarFieldNames_.size(); i++)
        {
            volScalarFieldList_[i].set
            (
	            timeIndexList_(recTime.timeName()),
                new volScalarField
                (
                    IOobject
                    (
                        volScalarFieldNames_[i],
                        recTime.timePath(),
                        base_.mesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    base_.mesh()
                )
	        );
	    }
	    
	    for(int i=0; i<volVectorFieldNames_.size(); i++)
	    {
            volVectorFieldList_[i].set
            (
	            timeIndexList_(recTime.timeName()),
                new volVectorField
                (
                    IOobject
                    (
                        volVectorFieldNames_[i],
                        recTime.timePath(),
                        base_.mesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    base_.mesh()
                )
	        );
	    }
	    
	    for(int i=0; i<surfaceScalarFieldNames_.size(); i++)
	    {
            surfaceScalarFieldList_[i].set
            (
                timeIndexList_(recTime.timeName()),
                new surfaceScalarField
                (
                    IOobject
                    (
                        surfaceScalarFieldNames_[i],
                        recTime.timePath(),
                        base_.mesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    base_.mesh()
                )
            );
        } 
    }
    Info << "Reading fields done" << endl;
}


void standardRecModel::readTimeSeries()
{
    bool firsttime = true;
    // fill the data structure for the time indices
    for (instantList::iterator it=timeDirs.begin(); it != timeDirs.end(); ++it)
    {
    	// set run-time
    	recTime.setTime(*it, it->value());
    	
    	
    	// skip constant
    	if (recTime.timeName() == "constant")
    	{
        	continue;
        }
        
        // skip zero
        if (skipZero_ and recTime.timeName() == "0")
        {
            if (verbose_)
    	    {
    	        Info << " ... skipping 0 in readTimeSeries()" << endl;
    	    }
            
            continue;
        }
        
        if (firsttime)
        {
            firsttime = false;
            recStartTime_ = recTime.value();
        }
	    recEndTime_ = recTime.value();
        
        // insert the time name into the hash-table with a continuous second index
        timeIndexList_.insert(recTime.timeName(), contTimeIndex);
        
        
        if (verbose_)
    	{
    		Info << "current time " << recTime.timeName() << endl;
    		Info << "insert " << recTime.timeName() << " , " << contTimeIndex << endl;
    	}
        
        
        // insert the time value
        timeValueList_.insert(contTimeIndex, recTime.timeOutputValue());
        
        // increment continuousIndex
        contTimeIndex++;
        
        if (verbose_)
    	{
        	Info << "contTimeIndex " << contTimeIndex << endl;
        }
    }

    if (verbose_)
    {
        Info << endl;
        Info << "Found " << label(timeDirs.size()) << " time folders" << endl;
        Info << "Found " << label(timeIndexList_.size()) << " time steps" << endl;
        Info << "database start time = " << recStartTime_ << endl;
        Info << "database end time = " << recEndTime_ << endl;
    }
}


void standardRecModel::exportVolScalarField(word fieldname, volScalarField& field)
{
    field = exportVolScalarField(fieldname, virtualTimeIndex);
}


void standardRecModel::exportVolVectorField(word fieldname, volVectorField& field)
{
    field = exportVolVectorField(fieldname, virtualTimeIndex);
}


void standardRecModel::exportSurfaceScalarField(word fieldname, surfaceScalarField& field)
{
    field = exportSurfaceScalarField(fieldname, virtualTimeIndex);
}


const volScalarField& standardRecModel::exportVolScalarField(word fieldname, label index)
{
    const label fieldI = getVolScalarFieldIndex(fieldname, index);
    
    return volScalarFieldList_[fieldI][index];
}

const volVectorField& standardRecModel::exportVolVectorField(word fieldname, label index)
{
    const label fieldI = getVolVectorFieldIndex(fieldname, index);
    
    return volVectorFieldList_[fieldI][index];
}

const surfaceScalarField& standardRecModel::exportSurfaceScalarField(word fieldname, label index)
{
    const label fieldI = getSurfaceScalarFieldIndex(fieldname, index);
    
    return surfaceScalarFieldList_[fieldI][index];
}


SymmetricSquareMatrix<scalar>& standardRecModel::recurrenceMatrix()
{
    return recurrenceMatrix_;
}


const HashTable<label,word>& standardRecModel::timeIndexList() const
{
    return timeIndexList_;
}


label standardRecModel::lowerSeqLim() const
{
    return lowerSeqLim_; 
}

label standardRecModel::upperSeqLim() const
{
    return upperSeqLim_; 
}


label standardRecModel::numRecFields() const
{
    return numRecFields_;
}

label standardRecModel::numDataBaseFields() const
{
    return numRecFields_;
}


void standardRecModel::updateRecFields()
{
    virtualTimeIndex = virtualTimeIndexNext;
    virtualTimeIndexNext++;
    
    if (virtualTimeIndexNext > sequenceEnd)
    {
        virtualTimeIndexListPos++; // this is problematic with noPath
        
        sequenceStart = virtualTimeIndexList_[virtualTimeIndexListPos].first();
        sequenceEnd = virtualTimeIndexList_[virtualTimeIndexListPos].second();
        
        virtualTimeIndexNext = sequenceStart;
        
        if (verbose_)
        {
            Info << " new sequence (start/end) : " << sequenceStart << " / " << sequenceEnd << endl;
        }
    }

    if (verbose_)
    {
        Info << "\nUpdating virtual time index to " << virtualTimeIndex << ".\n" << endl;
    }
}


void standardRecModel::writeRecMatrix() const
{
    OFstream matrixFile("recurrenceMatrix");
    matrixFile << recurrenceMatrix_;
}








// tmp<surfaceScalarField> standardRecModel::exportAveragedSurfaceScalarField(word fieldname, scalar threshold, label index)
// {
//     label timeIndex;
//     if (index < 0)
//     {
//         timeIndex = virtualTimeIndex;
//     }
//     else
//     {
//         timeIndex = index; 
//     }
//     const label fieldI = getSurfaceScalarFieldIndex(fieldname, timeIndex);
//     
//     tmp<surfaceScalarField> tAveragedSurfaceScalarField(surfaceScalarFieldList_[fieldI][timeIndex]);
//     
//     label counter = 1;
//     scalar recErr;
//     label delay = 10;
//     label lastMin = -1000;
//     
//     for(int runningTimeIndex = 1; runningTimeIndex < numRecFields_-1 ; runningTimeIndex++)
//     {
//         recErr = recurrenceMatrix_[timeIndex][runningTimeIndex];
//         if(recErr > threshold) continue;
//         if(recErr > recurrenceMatrix_[timeIndex][runningTimeIndex-1]) continue;
//         if(recErr > recurrenceMatrix_[timeIndex][runningTimeIndex+1]) continue;
//         if(abs(runningTimeIndex - timeIndex) < delay) continue;
//         if(abs(runningTimeIndex - lastMin) < delay) continue;
// 
//         lastMin = runningTimeIndex;
//         counter++;
//         tAveragedSurfaceScalarField += surfaceScalarFieldList_[fieldI][runningTimeIndex];
//     }
//     
//     tAveragedSurfaceScalarField /= counter;
//     return tAveragedSurfaceScalarField;
// }

void standardRecModel::exportAveragedVolVectorField(volVectorField& smoothfield, word fieldname, scalar threshold, label index) const
{
    label timeIndex;
    if (index < 0)
    {
        timeIndex = virtualTimeIndex;
    }
    else
    {
        timeIndex = index; 
    }
    const label fieldI = getVolVectorFieldIndex(fieldname, timeIndex);
    
    smoothfield = volVectorFieldList_[fieldI][timeIndex];
     
    label counter = 1;
    scalar recErr;
    label delay = 1;
    label lastMin = -1000;

    for(int runningTimeIndex = 0; runningTimeIndex < numRecFields_ ; runningTimeIndex++)
    {
        recErr = recurrenceMatrix_[timeIndex][runningTimeIndex];
        if(recErr > threshold) continue;
   //     if(recErr > recurrenceMatrix_[timeIndex][runningTimeIndex-1]) continue;
   //     if(recErr > recurrenceMatrix_[timeIndex][runningTimeIndex+1]) continue;
        if(abs(runningTimeIndex - timeIndex) < delay) continue;
        if(abs(runningTimeIndex - lastMin) < delay) continue;

        lastMin = runningTimeIndex;
        counter++;
        smoothfield += volVectorFieldList_[fieldI][runningTimeIndex];
    }
    Info << "time index = " << index << ", counter = " << counter << endl;
    smoothfield /= counter;
}








// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
