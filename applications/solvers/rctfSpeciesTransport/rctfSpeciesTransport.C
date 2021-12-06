/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling

    Contributing authors:
    Thomas Lichtenegger, Gerhard Holzinger, Sanaz Abbasi
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
    Turbulent Transport Recurrence Solver for modal decomposition

Description
    Solves a transport equation for a passive scalar on a single-phase solution
    for a solver based on recurrence statistics

Rules
    Solution data to compute the recurrence statistics from, needs to
    reside in $CASE_ROOT/dataBase(0...N)
    Time step data in the first dataBase needs to be evenly spaced in time
    A list of indices for the corresponding incoherent fields to coherent ones
    should be provided.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "fvOptions.H"

#include "recBase.H"
#include "recModel.H"

#include "clockModel.H"

#include "objectRegistry.H"
#include "VectorSpace.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    scalar relaxCoeff(0.0);

    //create recBases according to a list of recProperties
    #include "createRecBase.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    label recTimeIndex(0);
    label currTimeIndex(0);

    scalar recTimeStep_=recBases[0].recM().recTimeStep();
    labelPairList incPairTimeIndex_(0);

    IFstream pairFile("incIndexPairList");
    pairFile >> incPairTimeIndex_;

    while (runTime.run())
    {

        myClock().start(1,"Global");
        runTime++;

        myClock().start(11,"Total");

        Info<< "Time = " << runTime.timeName() << nl << endl;

        myClock().start(2,"fieldUpdate");

        if ( runTime.timeOutputValue() - (recTimeIndex+1)*recTimeStep_ + 1.0e-5 > 0.0 )
        {
            Info<< "Updating fields at run time " << runTime.timeOutputValue()
                << " corresponding to recurrence time " << (recTimeIndex+1)*recTimeStep_ << ".\n" << endl;
            recBases[0].updateRecFields();
            #include "readFields.H"

            recTimeIndex++;
        }

        myClock().stop("fieldUpdate");

        #include "continuityErrCalc.H"

        myClock().start(3,"speciesEqn");
        #include "CEq.H"
        myClock().stop("speciesEqn");

        myClock().stop("Total");

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        myClock().stop("Global");

    }


    myClock().evalPar();
    myClock().normHist();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
