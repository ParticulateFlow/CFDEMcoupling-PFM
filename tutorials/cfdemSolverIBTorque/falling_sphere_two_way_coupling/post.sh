#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine
# run settlingTest CFD part
# Christoph Goniva - Feb. 2011
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_CFDEMIB_fixed_RS"
logfileName="log_$headerText"
solverName="cfdemSolverIB"
nrProcs="4"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
postproc="true"
#--------------------------------------------------------------------------------#

#- get VTK data from liggghts dump file
cd $casePath/DEM/post
lpp dump*

#- get VTK data from CFD sim
cd $casePath/CFD
reconstructParMesh -constant -mergeTol 1e-06
reconstructPar -noLagrangian
foamToVTK

#- start paraview
paraview

