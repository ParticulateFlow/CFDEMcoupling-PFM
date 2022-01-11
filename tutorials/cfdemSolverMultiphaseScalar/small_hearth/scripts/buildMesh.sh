#!/bin/bash

#===================================================================#
# construct mesh
#===================================================================#

#--------------------------------------------------------------------------------#
#- run settings file
scriptPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
. $scriptPath/scriptSettings.sh $1
#--------------------------------------------------------------------------------#

cd $casePath/CFD

# blockmesh
if [ -f "$casePath/CFD/constant/polyMesh/blockMeshDict" ]; then
    echo "Meshing with blockMesh"
    blockMesh
fi
# snappyhexmesh
if [ -f "$casePath/CFD/system/snappyHexMeshDict" ]; then
    echo "Meshing with snappyHexMesh"
    surfaceFeatureExtract
    decomposePar
    mpirun -np $nrProcs snappyHexMesh -overwrite -parallel
    reconstructParMesh -constant -mergeTol 1e-6
    rm -r processor*
fi

# UNV to foam
if ls $casePath/CFD/constant/polyMesh/*.unv 1> /dev/null 2>&1; then
    echo "Meshing with ideasUnvToFoam"
    ideasUnvToFoam $casePath/CFD/constant/polyMesh/*.unv
fi

cd $currentPath

	




