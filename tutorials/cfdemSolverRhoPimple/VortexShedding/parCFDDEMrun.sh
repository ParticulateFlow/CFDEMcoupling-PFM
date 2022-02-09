#!/bin/bash
#------------------------------------------------------------------------------
# parCFDDEMrun script for VortexShedding test case
# run VortexShedding CFD-DEM
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
headerText="run_parallel_cfdemSolverRhoPimple_VortexShedding"
logfileName="log_$headerText"
solverName="cfdemSolverRhoPimple"
nrProcs="4"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
#------------------------------------------------------------------------------

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

