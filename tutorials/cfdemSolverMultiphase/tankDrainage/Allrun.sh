#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine
# run tankDrainage
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
postProcessing=true

# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi

cd $casePath/CFD
cp -r 0.org 0
setFields

if [ -f "$casePath/DEM/post/restart/liggghts.restart" ];  then
    echo "LIGGGHTS init was run before - using existing restart file"
else
    #- run DEM in new terminal
    $casePath/parDEMrun.sh
fi

#- run parallel CFD-DEM in new terminal
bash $casePath/parCFDDEMrun.sh

if [ "$postProcessing" = true ]; then
    cd $casePath
    bash $casePath/postrun.sh
    octave $casePath/volFlow.m
fi
