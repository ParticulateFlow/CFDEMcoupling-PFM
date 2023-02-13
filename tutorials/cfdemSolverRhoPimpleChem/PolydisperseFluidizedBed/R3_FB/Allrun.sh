#!/bin/bash
#------------------------------------------------------------------------------
# Allrun script for fluidized bed R3 chemistry test case
# run R3_FB test case
# Daniel Queteschiner - March 2022
#------------------------------------------------------------------------------

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

#- check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi

if [ -f "$casePath/DEM/post/restart/liggghts.restart" ];  then
    echo "LIGGGHTS init was run before - using existing restart file"
else
    #- run parallel DEM
    $casePath/parDEMrun.sh
fi

#- run parallel CFD-DEM
bash $casePath/parCFDDEMrun.sh
