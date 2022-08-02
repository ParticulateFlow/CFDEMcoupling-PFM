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
    rcfdemSolverHeattransfer

Description
    Solves heat transfer between fluid and particles based on rCFD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "cfdemCloudRec.H"
#include "recBase.H"
#include "recModel.H"
#include "recPath.H"

#include "cfdemCloudEnergy.H"
#include "clockModel.H"
#include "thermCondModel.H"
#include "energyModel.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"

    cfdemCloudRec<cfdemCloudEnergy> particleCloud(mesh);
    recBase recurrenceBase(mesh);
    #include "updateFields.H"
    #include "updateRho.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nCalculating particle trajectories based on recurrence statistics\n" << endl;

    label recTimeIndex = 0;
    label stepCounter = 0;
    label recTimeStep2CFDTimeStep = recurrenceBase.recM().recTimeStep2CFDTimeStep();

    // control coupling behavior in case of substepping
    // assumes constant timestep size
    label counter = 0;
    label couplingSubStep = recurrenceBase.couplingSubStep();
    double dtProp =  particleCloud.dataExchangeM().couplingTime() / runTime.deltaTValue();
    label dtDEM2dtCFD = int(dtProp + 0.5);
    Info << "deltaT_DEM / deltaT_CFD = " << dtDEM2dtCFD << endl;
    if (dtDEM2dtCFD > 1)
        Info << "coupling at substep " << couplingSubStep << endl;


    while (runTime.run())
    {
        runTime++;

        // do stuff (every lagrangian time step)
        particleCloud.clockM().start(1,"Global");

        Info << "Time = " << runTime.timeName() << nl << endl;

        particleCloud.clockM().start(2,"Coupling");

        particleCloud.evolve(voidfraction,Us,URec);

        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");
        #include "updateRho.H"
        #include "TEqn.H"
        particleCloud.clockM().stop("Flow");

        stepCounter++;

        particleCloud.clockM().start(32,"ReadFields");
        if (stepCounter == recTimeStep2CFDTimeStep)
        {
            recurrenceBase.updateRecFields();
            #include "updateFields.H"
            recTimeIndex++;
            stepCounter = 0;
            recTimeStep2CFDTimeStep = recurrenceBase.recM().recTimeStep2CFDTimeStep();
        }
        particleCloud.clockM().stop("ReadFields");

        particleCloud.clockM().start(33,"Output");
        runTime.write();
        particleCloud.clockM().stop("Output");

        particleCloud.clockM().stop("Global");

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
