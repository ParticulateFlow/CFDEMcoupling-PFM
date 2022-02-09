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
headerText="run_parallel_cfdemSolverPiso_ErgunTestMPI_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverIB"
nrProcs="4"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
reconstructCase="false"  # true | false
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="true"
postproc="false"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

#- case needs special reconstruction
if [ $reconstructCase == "true" ]
  then
    cd $casePath/CFD
    reconstructParMesh -mergeTol 1e-06
    reconstructPar -noLagrangian
fi

if [ $runOctave == "true" ]
  then

    cd $casePath/CFD/octave
    octave postproc.m
    evince pos_z_two_part_rec_glow.eps
    evince vel_z_two_part_rec_glow.eps
    #display pos_y_two_part_rec_glow.png &
    #display vel_y_two_part_rec_glow.png &
fi

if [ $postproc == "true" ]
  then
    #- get VTK data from liggghts dump file
    cd $casePath/DEM/post
    python -i $CFDEM_LPP_DIR/lpp.py  dump.liggghts_run

    #- get VTK data from CFD sim
    cd $casePath/CFD
    foamToVTK

    #- start paraview
    paraview

    #- keep terminal open (if started in new terminal)
    echo "...press enter to clean up case"
    echo "press Ctr+C to keep data"
    read
fi

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

