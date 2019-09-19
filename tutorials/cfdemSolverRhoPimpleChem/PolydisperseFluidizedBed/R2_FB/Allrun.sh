#!/bin/bash

#===================================================================#
# allrun script for testcase
# M. Efe Kinaci - Sep 2018
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
export casePath

cd $casePath/CFD/
blockMesh

echo "Run Simulation"
bash $casePath/parCFDDEMrun.sh   
