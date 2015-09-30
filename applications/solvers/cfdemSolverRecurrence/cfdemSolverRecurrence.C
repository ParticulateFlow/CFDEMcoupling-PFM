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
#include "fvIOoptionList.H"	// no solver = no fvOptions

#include "mathematicalConstants.H" // provides: e, pi, 2pi and pi/2

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
    #include "createRecurrentFields.H"
    #include "createFields.H"
    #include "createFvOptions.H"	// no solver = no fvOptions
  
    cfdemCloudRec particleCloud(mesh);
    #include "checkModelType.H"

    //simpleControl simple(mesh);	// no solver = no solution control
    
    // random recurrence stuff
    #include "Random.H"	// random numbers
    Random ranGen(osRandomInteger());
    
    label sequenceLength(0);
    
    // minimum sequence length = the larger of 1 or a twentieth of all
    const label lowerSeqLim(max(1, label(timeIndexList.size()/20)));
    
    // maximum sequence length = a fifth of all
    const label upperSeqLim(label(timeIndexList.size()/5));
    
    sequenceLength = ranGen.integer(lowerSeqLim, upperSeqLim);
    
    scalar nextBestMinimum(GREAT);
    scalar secondBestMinimum(GREAT);
    label sequenceStart(0);
    label sequenceStart2(0);
    label sequenceStartOld(0);
    label virtualTimeIndex(0);
    label virtualStartIndex(0);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating particle trajectories based on recurrence statistics\n" << endl;
    
    label recTimeIndex(0);
    
    particleCloud.initRecFields();
    while (runTime.run())
    {
        // advance time
        runTime++;
        
        // do stuff (every lagrangian time step)
        particleCloud.clockM().start(1,"Global");

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // do particle stuff
        particleCloud.clockM().start(2,"Coupling");

        particleCloud.evolve(voidfraction,Us,U);


        #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");
        
        if ( runTime.timeOutputValue() / ((recTimeIndex+1)*dtCur) >= 1.0 )
        {
        	
        	// advance!
	  	particleCloud.updateRecFields();
		
        	virtualTimeIndex++;
        	recTimeIndex++;
        	
        	alpha2 = alpha2pl[virtualTimeIndex];
        	U2 = U2pl[virtualTimeIndex];
        	
        	if (!verbose)
        	{
        		Info << nl << "Time = " << runTime.timeName() << endl;
        		Info<< "ExecutionTime = "
            		<< runTime.elapsedCpuTime()
            		<< " s\n\n" << endl;
        	}
        	
        	
        	timeSequenceFile << virtualTimeIndex << endl;
        }
        
        
        
        #include "recurrenceTreatment.H"
        
      

       
        
        
        // write stuff at output time
        if (runTime.outputTime())
        {        	
        //	magicParticles.writeConcField();
        	
        	alpha2.write();
        	U2.write();
        }
        
        runTime.write(); // does not work, yet
	
	

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
        
        if (verbose)
        {
        	 Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        }
    }
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
