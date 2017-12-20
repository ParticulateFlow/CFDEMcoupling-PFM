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

Application
    rSmoother

Description
    Loops over all recurrence times and averages fields over given similarity range


\*---------------------------------------------------------------------------*/

 #include "fvCFD.H"
 #include "fvOptions.H"

#include "recBase.H"
#include "recModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
  
 
    recBase recurrenceBase(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // set time step to that of recurrence database
    runTime.setDeltaT(recurrenceBase.recM().recTimeStep());

    // check start and end time
    if (abs(runTime.startTime().value() - recurrenceBase.recM().recStartTime()) > 1e-5)
    {
        Info << "Stopping. Start time and database start time are different." << endl;
	Info << "Start time = " << runTime.startTime().value() << endl;
	Info << "Database start time = " << recurrenceBase.recM().recStartTime() << endl;
        return 0;
    }

    if (runTime.endTime().value() > recurrenceBase.recM().recEndTime())
    {
        runTime.setEndTime(recurrenceBase.recM().recEndTime());
        Info << "End time set to database end time." << endl;
    }

    label index = -1;

    Info<< "\nSmoothing recurrence statistics\n" << endl;
    
    while (runTime.run())
    {
        // runtime can't be larger than recurrence database size
        index = runTime.timeIndex();
	#include "updateFields.H"
        runTime++;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
