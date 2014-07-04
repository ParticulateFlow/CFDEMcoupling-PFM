#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run packedBedTemp
# Christoph Goniva - June 2014
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi

# check if DEM case was run
if [ -f "$casePath/DEM/liggghts.restart" ]; then
    echo "DEM restart file found"
else
    echo "starting DEM run..."
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$casePath"
    headerText="run_liggghts_packedBedTemp_DEM"
    logfileName="log_$headerText"
    solverName="in.liggghts_init"
    nrProcs=4
    machineFileName="none"
    debugMode="off"
    #--------------------------------------------------------------------------------#

    #- clean up case
    rm -r $casePath/DEM/post/*

    #- call function to run DEM case
    parDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode
fi

#- run parallel CFD-DEM in new terminal
gnome-terminal --title='cfdemSolverPisoScalar packedBedTemp CFD'  -e "bash $casePath/parCFDDEMrun.sh" 

