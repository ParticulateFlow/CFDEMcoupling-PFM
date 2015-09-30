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
#include "recModel.H"
#include <unistd.h>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(recModel, 0);

defineRunTimeSelectionTable(recModel, dictionary);


// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::recModel::recModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    recTime("dataBase", "", "../system", "../constant", false),
    timeDirs(recTime.times()),
    timeIndexList(label(timeDirs.size())-1),
    timeValueList(label(timeDirs.size())-1),
    contTimeIndex(0),
    sequenceStart(0),
    sequenceStart2(0),
    sequenceStartOld(0),
    virtualTimeIndex(0),
    virtualStartIndex(0)
{
    if (verbose)
    {
    	// be informative on properties of the "recTime" Time-object
    	Info << "recTime.rootPath() " << recTime.rootPath() << endl;
    	Info << "recTime.caseName() " << recTime.caseName() << endl;
    	Info << "recTime.path() " << recTime.path() << endl;
    	Info << "recTime.timePath() " << recTime.timePath() << endl;
	Info << "recTime.timeName() " << recTime.timeName() << endl;
	Info << "timeDirs " << timeDirs << endl;
    }
    readTimeSeries();
    checkTimeStep();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::recModel::~recModel()
{}

// * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * //

virtual void Foam::recModel::updateRecFields()
{}

// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * //

void Foam::recModel::readTimeSeries()
{
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
        
        
        // insert the time name into the hash-table with a continuous second index
        timeIndexList.insert(recTime.timeName(), contTimeIndex);
        
        
        if (verbose)
    	{
    		Info << "current time " << recTime.timeName() << endl;
    		Info << "insert " << recTime.timeName() << " , " << contTimeIndex << endl;
    	}
        
        
        // insert the time value
        timeValueList.insert(contTimeIndex, recTime.timeOutputValue());
        
        // increment continuousIndex
        contTimeIndex++;
        
        if (verbose)
    	{
        	Info << "contTimeIndex " << contTimeIndex << endl;
        }
    }

    if (verbose)
    {
    	Info << endl;
    	Info << "Found " << label(timeDirs.size()) << " time folders" << endl;
    	Info << "Used " << (label(timeDirs.size())-1) << " time folders" << endl;
        Info << "Found " << label(timeIndexList.size()) << " time steps" << endl;
    }
}

void Foam::recModel::checkTimeStep()
{
   // check time step of provided data
    scalar dtCur(0.0);
    scalar dtOld(0.0);
    
    if (verbose)
    {
    	Info << "timeValueList : " << timeValueList << endl;
    }
    
    forAll(timeValueList, i)
    {
    	// compute time step
    	if (timeDirs[i].value() == timeDirs.last().value())
    	{
    		if (verbose)
    		{
    			Info << ".. leaving loop at " << timeDirs[i] << endl;
    		}
    		// leave loop
    		break;
    	}
    	
    	if (verbose)
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
	
	if (verbose)
    {
		Info << "Setting deltaT to " << dtCur << endl;
		Info << "Actual recTime.deltaT = " << recTime.deltaTValue() << endl;
		Info << "Actual runTime.deltaT = " << runTime.deltaTValue() << endl;
    }
}

void Foam::recModel::writeRecMatrix()
{
    // create output file
    std::ostringstream str_pid;
    str_pid << pid();
    const std::string my_pid(str_pid.str());
    OFstream matrixFile("recurrenceMatrix."+my_pid);
     // write recurrence matrix
    matrixFile << recurrenceMatrix;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
