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

Application
    cfdemSolverRecurrence

Description
    Solves a transport equation for a passive scalar on a two-phase solution
    Test-bed for a solver based on recurrence statistics

Rules
	Solution data to compute the recurrence statistics from, needs to 
		reside in $CASE_ROOT/dataBase
	Time step data in dataBase needs to be evenly spaced in time

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvIOoptionList.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "cfdemCloudRec.H"
#include "clockModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    cfdemCloudRec particleCloud(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating particle trajectories based on recurrence statistics\n" << endl;
    
    label recTimeIndex(0);
    scalar recTimeStep_=particleCloud.recM().recTimeStep();
    
    while (runTime.run())
    {
        runTime++;
        
        // do stuff (every lagrangian time step)
        particleCloud.clockM().start(1,"Global");

        Info<< "Time = " << runTime.timeName() << nl << endl;
	
        #include "TEq.H"

        particleCloud.clockM().start(26,"Flow");
        
        if ( runTime.timeOutputValue() - (recTimeIndex+1)*recTimeStep_ + 1.0e-5 > 0.0 )
        {
	    Info << "Updating fields at run time " << runTime.timeOutputValue()
	        << " corresponding to recurrence time " << (recTimeIndex+1)*recTimeStep_ << ".\n" << endl;
            particleCloud.updateRecFields();
	    #include "readFields.H"
            recTimeIndex++;
        }

        runTime.write();	

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
        
    }
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
