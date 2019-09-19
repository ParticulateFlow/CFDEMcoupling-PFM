#!/bin/bash

#===================================================================#
# allrun script for testcase
# M. Efe Kinaci - Sep 2018
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
export casePath

cd $casePath/CFD/
blockMesh

if [ -f "$casePath/DEM/post/restart/liggghts.restart" ]; then
    echo "LIGGGHTS init was run before - using existing restart file"
else
    #- run DEM in new terminal
    $casePath/parDEMrun.sh
fi

echo "Run Simulation"
bash $casePath/parCFDDEMrun.sh
