#!/bin/bash
#------------------------------------------------------------------------------
# Allrun script for HeatTransferRanzMarshall test case
# run HeatTransferRanzMarshall
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
fi

#- run parallel CFD-DEM in new terminal
bash $casePath/parCFDDEMrun.sh
