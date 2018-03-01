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

# check if initialization was done
if [ -f "$casePath/DEM/post/restart/liggghts.restart" ];  then
    echo "Initialization was run before - using existing restart file"
else
	cd $casePath
	./parDEMrun.sh
fi

cd $casePath/CFD
cp -r 0.org 0
setFields

cd $casePath
./parCFDDEMrun.sh

# generate files for post processing
cd $casePath
./postRun.sh
