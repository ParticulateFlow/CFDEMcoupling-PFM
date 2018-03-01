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
if [ -f "$casePath/DEM/post/restart/liggghts.restartCFDEM_5.000000" ];  then
    echo "Initialization was run before - using existing restart file"
else
	cd $casePath
	bash init.sh
fi

# run the drainage simulation
cd $casePath/CFD
cp system/controlDict.run system/controlDict
cp constant/couplingProperties.run constant couplingProperties
cp -r $casePath/DEM/post/restart/5/ .
decomposePar -force
mpirun -np $nrProcs cfdemSolverMultiphase -parallel |& tee $casePath/log_run

# generate files for post processing
cd $casePath
./postRun.sh
