#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine
# run settlingTest CFD part
# Christoph Goniva - Feb. 2011
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_cfdemSolverPisoMS_ErgunTestMPI_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPimpleDyMMS_22x"
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"

cleanUp="true"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

#------------------------------#
#  octave

#- change path
cd octave

#- rmove old graph
rm cfdemSolverPisoMS_ErgunTestMPI.eps

#- run octave
octave totalPressureDrop.m

#- show plot
evince cfdemSolverPisoMS_ErgunTestMPI.eps
#------------------------------#

#- copy log file to test harness
cp ../../$logfileName $testHarnessPath
cp cfdemSolverPisoMS_ErgunTestMPI.eps $testHarnessPath

if [ $cleanUp == "true" ]
  then
    #- clean up case
    echo "deleting data at: $casePath :\n"
    source $WM_PROJECT_DIR/bin/tools/CleanFunctions
    cd $casePath/CFD
    cleanCase
    rm -r $casePath/CFD/clockData
    rm -r $casePath/DEM/post/*
    (cd $casePath/DEM/post && touch dummy)
    echo "done"
fi

