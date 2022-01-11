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
    rcfdemSolverForcedTracers

Description
    Moves tracers according to the activated force models on pressure and velocity
    fields provided by a recurrence process

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "fvOptions.H"

#include "recBase.H"
#include "recModel.H"

#include "cfdemCloud.H"
#include "clockModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"

    cfdemCloud particleCloud(mesh);
    recBase recurrenceBase(mesh);

    const IOdictionary& recProps = mesh.lookupObject<IOdictionary>("recProperties");
    bool useRecP(recProps.lookupOrDefault<bool>("useRecP",false));
    bool useRecK(recProps.lookupOrDefault<bool>("useRecK",false));

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nCalculating particle trajectories based on recurrence statistics\n" << endl;

    label recTimeIndex = 0;
    label stepCounter = 0;
    label recTimeStep2CFDTimeStep = recurrenceBase.recM().recTimeStep2CFDTimeStep();

    while (runTime.run())
    {
        runTime++;

        // do stuff (every lagrangian time step)
        particleCloud.clockM().start(1,"Global");

        Info << "Time = " << runTime.timeName() << nl << endl;

        particleCloud.clockM().start(2,"Coupling");

        particleCloud.evolve(voidfraction,Us,URec);

        particleCloud.clockM().stop("Coupling");

        stepCounter++;

        if (stepCounter == recTimeStep2CFDTimeStep)
        {
            recurrenceBase.updateRecFields();
            #include "updateFields.H"
            recTimeIndex++;
            stepCounter = 0;
            recTimeStep2CFDTimeStep = recurrenceBase.recM().recTimeStep2CFDTimeStep();
        }

        particleCloud.clockM().start(27,"Output");
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
