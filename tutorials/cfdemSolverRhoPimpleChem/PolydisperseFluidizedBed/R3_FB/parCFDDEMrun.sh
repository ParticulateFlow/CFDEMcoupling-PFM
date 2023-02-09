#!/bin/bash
#------------------------------------------------------------------------------
# parCFDDEMrun script for R3_FB test case
# run R3_FB CFD-DEM
# Daniel Queteschiner - March 2022
#------------------------------------------------------------------------------

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#------------------------------------------------------------------------------
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="R3_FluidizedBed"
logfileName="log_$headerText"
solverName="cfdemSolverRhoPimpleChem"
nrProcs="4"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="true"

#------------------------------------------------------------------------------

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode


if [ $runOctave == "true" ]
    then
        #------------------------------#
        #  octave

        #- change path
        cd octave

        #- run octave
        octave plotData.m
fi

