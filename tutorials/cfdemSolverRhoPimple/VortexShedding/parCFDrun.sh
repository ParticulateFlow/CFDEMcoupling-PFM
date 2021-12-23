#!/bin/bash
#------------------------------------------------------------------------------
# parCFDrun script for VortexShedding test case
# run VortexShedding CFD
# Daniel Queteschiner - November 2021
#------------------------------------------------------------------------------

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#------------------------------------------------------------------------------
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_rhoPimpleFoam_VortexShedding"
logfileName="log_$headerText"
solverName="rhoPimpleFoam"
nrProcs="4"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
#------------------------------------------------------------------------------

#- call function to run a parallel CFD-DEM case
parCFDrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

