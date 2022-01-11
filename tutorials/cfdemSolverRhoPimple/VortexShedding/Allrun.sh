#!/bin/bash
#------------------------------------------------------------------------------
# Allrun script for VortexShedding test case
# run VortexShedding
# Daniel Queteschiner - November 2021
#------------------------------------------------------------------------------

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

#- check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
    transformPoints -scale "(0.001 0.001 0.001)"
    checkMesh
    #- make the linear system more diagonal dominant to speed-up the linear solvers
    renumberMesh -overwrite -noFunctionObjects
fi

#- run parallel CFD
#bash $casePath/parCFDrun.sh
#- run parallel CFD-DEM
bash $casePath/parCFDDEMrun.sh
