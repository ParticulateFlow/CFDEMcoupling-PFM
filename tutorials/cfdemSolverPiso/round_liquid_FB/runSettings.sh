#!/bin/bash

#===================================================================#
# run settings for case
# Tim MJ Nijssen - May 2018
#===================================================================#

#- decide what run settings to use
if [ -z "$1" ]; then
    #- ask runtype if not specified
    echo "Please specify runtype (DEM/CFDDEM):"
    read runType
else
    #- use specified type
    runType=$1
fi

#- check runType
if [ "$runType" == DEM ] || [ "$runType" == CFDDEM ]; then
    echo "Running settings for $runType"
else
    echo "Error: invalid runType"
    read
fi

if [ "$runType" == DEM ]; then
    #- DEM run settings
    nrProcs="16"
    runTitle="DEMrun"
    headerText="$runTitle"
    logfileName="log_DEM"
    solverName="in.liggghts_init"
    machineFileName="none"   # yourMachinefileName | none
    debugMode="off"
fi

if [ "$runType" == CFDDEM ]; then
    #- CFDDEM run settings
    nrProcs="16"
    headerText="$runTitle"
    logfileName="log_CFDDEM"
    solverName="cfdemSolverPiso"
    machineFileName="none"   # yourMachinefileName | none
    debugMode="off"          # on | off| strict
    testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
fi
