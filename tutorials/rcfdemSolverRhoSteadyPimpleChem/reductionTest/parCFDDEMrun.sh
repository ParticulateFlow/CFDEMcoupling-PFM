#!/bin/bash

#===================================================================#
# CFD-DEM run script for ore reduction testcase
# Thomas Lichtenegger - January 2023
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_rcfdemSolverRhoSteadyPimpleChem_reductionTest_CFDDEM"
logfileName="log_$headerText"
solverName="rcfdemSolverRhoSteadyPimpleChem"
nrProcs="8"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode
