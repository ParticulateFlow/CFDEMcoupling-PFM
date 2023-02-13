#!/bin/bash
#------------------------------------------------------------------------------
# allrun script for testcase
# M. Efe Kinaci - Sep 2018
#------------------------------------------------------------------------------

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#------------------------------------------------------------------------------
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="R2_FluidBed"
logfileName="log_$headerText"
solverName="cfdemSolverRhoPimpleChem"
nrProcs="4"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="true"
postproc="false"

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
    
if [ $postproc == "true" ]
  then

    #- keep terminal open (if started in new terminal)
    echo "simulation finished? ...press enter to proceed"
    read

    #- get VTK data from liggghts dump file
    cd $casePath/DEM/post
    python -i $CFDEM_LPP_DIR/lpp.py dump*.liggghts_run

    #- get VTK data from CFD sim
    cd $casePath/CFD
    reconstructPar

    #- start paraview
    paraview
 
    #- keep terminal open (if started in new terminal)
    echo "...press enter to clean up case"
    echo "press Ctr+C to keep data"
    read
fi

