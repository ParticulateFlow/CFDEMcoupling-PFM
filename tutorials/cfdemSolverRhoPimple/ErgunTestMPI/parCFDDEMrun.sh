#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine
# run Ergun test CFD part
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_cfdemSolverRhoPimple_ErgunTestMPI_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverRhoPimple"
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
        #------------------------------#
        #  octave

        #- change path
        cd octave

        #- rmove old graph
        rm cfdemSolverRhoPimple_ErgunTestMPI.eps

        #- run octave
        octave totalPressureDrop.m

        #- show plot
        evince cfdemSolverRhoPimple_ErgunTestMPI.eps

        #- copy log file to test harness
        cp ../../$logfileName $testHarnessPath
        cp cfdemSolverRhoPimple_ErgunTestMPI.eps $testHarnessPath
fi

if [ $postproc == "true" ]
  then

    #- keep terminal open (if started in new terminal)
    echo "simulation finisehd? ...press enter to proceed"
    read

    #- get VTK data from liggghts dump file
    cd $casePath/DEM/post
    python -i $CFDEM_LPP_DIR/lpp.py dump*.liggghts_run

    #- get VTK data from CFD sim
    cd $casePath/CFD
    foamToVTK                                                   #- serial run of foamToVTK

    #- start paraview
    paraview

    #- keep terminal open (if started in new terminal)
    echo "...press enter to clean up case"
    echo "press Ctr+C to keep data"
    read

fi

#- clean up case
echo "deleting data at: $casePath :\n"
source $WM_PROJECT_DIR/bin/tools/CleanFunctions
cd $casePath/CFD
#cleanCase
rm -r $casePath/CFD/clockData
rm $casePath/DEM/post/*.*
touch $casePath/DEM/post/.gitignore
rm $casePath/DEM/post/restart/*.*
rm $casePath/DEM/post/restart/liggghts.restartCFDEM*
touch $casePath/DEM/post/restart/.gitignore
echo "done"
