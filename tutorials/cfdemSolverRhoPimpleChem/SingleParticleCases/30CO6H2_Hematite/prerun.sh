#!/bin/bash

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/boundary" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi