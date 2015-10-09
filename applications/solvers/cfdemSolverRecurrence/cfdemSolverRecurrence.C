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
	Hardcoded Eulerian set-up
		The species is part of/ transported by the continuous phase
		"water" is the continuous phase
		"air" is the dispersed
		Make a donation to the get selectable phase names implemented

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
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
    cfdemCloudRec particleCloud(mesh);
    #include "checkModelType.H"

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

      //  #include "CourantNo.H"

        particleCloud.clockM().start(2,"Coupling");

        particleCloud.evolve(voidfraction,Us,URec);

        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");
        
        if ( runTime.timeOutputValue()  >= (recTimeIndex+1)*recTimeStep_ )
        {
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
