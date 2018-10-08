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
headerText="run_parallel_CFDEMIB_two_way_coupling"
logfileName="log_$headerText"
solverName="cfdemSolverIBTorque"
nrProcs="4"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="true"
postproc="false"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

if [ $runOctave == "true" ]
  then

    cd $casePath/CFD/octave
    octave plot_data.m
    evince angular_velocity_compare.eps
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

#- copy log file to test harness
cp ../../$logfileName $testHarnessPath

#- clean up case
echo "deleting data at: $casePath"
source $WM_PROJECT_DIR/bin/tools/CleanFunctions
cd $casePath/CFD
cleanCase
rm -r $casePath/CFD/clockData
rm $casePath/DEM/post/*.*
#rm -r $casePath/DEM/post/restart/*.*
#- preserve post directory
touch $casePath/DEM/post/.gitignore
echo "done"

