#!/bin/bash

#===================================================================#
# allrun script for testcase
# M. Efe Kinaci - Sep 2019
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
export casePath

#cd $casePath/CFD
#blockMesh

#$casePath/parDEMrun.sh

#bash $casePath/parCFDDEMrun.sh

export casePath
# check if mesh was built
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
    #- run DEM in new terminal
    $casePath/parDEMrun.sh
fi

#echo "Run Simulation"
#cd $casePath/CFD
#decomposePar    
#mpirun -np $nrProcs $solverName -parallel | tee -a $logfileName 

#- run parallel CFD-DEM in new terminal
#gnome-terminal -e "bash $casePath/parCFDDEMrun.sh"
bash $casePath/parCFDDEMrun.sh   
