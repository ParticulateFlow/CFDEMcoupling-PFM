#!/bin/bash
#------------------------------------------------------------------------------
# parCFDDEMrun script for redBloodCellPoiseuilleFlow test case
# run redBloodCellPoiseuilleFlow CFD-DEM
# Achuth N. Balachandran Nair - October 2018
#------------------------------------------------------------------------------

#- source CFDEM env vars
. ~/.bashrc

#- source CFDEM functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#------------------------------------------------------------------------------
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_cfdemSolverIBContinuousForcing_redBloodCellPoiseuilleFlow"
logfileName="log_$headerText"
solverName="cfdemSolverIBContinuousForcing"
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runPython="false"
postproc="false"
#------------------------------------------------------------------------------

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

if [ $runPython == "true" ]
  then
    cd $casePath
    python results.py
fi

if [ $postproc == "true" ]
  then
    #- get VTK data from liggghts dump file
    cd $casePath/DEM/post
    lpp dump*

    #- get VTK data from CFD sim
    cd $casePath/CFD
    reconstructParMesh -constant -mergeTol 1e-06
    reconstructPar -noLagrangian
    foamToVTK

    #- start paraview
    paraview

    #- keep terminal open (if started in new terminal)
    echo "...press enter to clean up case"
    echo "press Ctr+C to keep data"
    read
fi


