#!/bin/bash

#===================================================================#
# This tutorial initializes a drainage simulation of a small-scale
# blast furnace hearth. A short coupled CFD-DEM simulation must be
# carried out so the particle bed can settle in buoyancy equilibrium.
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

##########################################################################################
# The following commands were used to set up the run case.
# The outlet boundary conditions were manually changed in ../run/CFD/1.org/p_rgh to 
# type fixedValue; value uniform 0; and in ../run/CFD/1.org/U to type zeroGradient;
##########################################################################################

#cd $casePath
#cp ./CFD/0.org/pSmooth ./CFD/1
#cp -r ./CFD/1/ ../run/CFD/1.org
#cp ./DEM/post/restart/liggghts.restartCFDEM_1.000000 ../run/DEM/post/restart/liggghts.restart

