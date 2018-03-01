#!/bin/bash

#===================================================================#
# allrun script for cfdemSolverMultiphase 
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
nrProcs=4;

# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
	cd $casePath/CFD
    ./mesh.sh
fi

#run the drainage simulation
cd $casePath
cp -r ./CFD/1.org ./CFD/1
./parCFDDEMrun.sh

# generate files for post processing
cd $casePath
./postRun.sh
