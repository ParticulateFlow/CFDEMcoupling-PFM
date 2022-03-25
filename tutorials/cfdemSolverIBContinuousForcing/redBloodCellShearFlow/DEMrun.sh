#!/bin/bash
#------------------------------------------------------------------------------
# DEMrun script for redBloodCellShearFlow test case
# init redBloodCellShearFlow
# Achuth N. Balachandran Nair - October 2018
#------------------------------------------------------------------------------

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#------------------------------------------------------------------------------
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath="$casePath"
headerText="run_liggghts_init_DEM"
logfileName="log_$headerText"
solverName="in.liggghts_init"
debugMode="off"
#------------------------------------------------------------------------------

#- call function to run DEM case
DEMrun $logpath $logfileName $casePath $headerText $solverName $debugMode

