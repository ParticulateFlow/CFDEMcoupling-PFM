#!/bin/bash

#===================================================================#
# DEMrun script for case (init)
# Tim MJ Nijssen - September 2021
# based on: Daniel Queteschiner - June 2014
#===================================================================#

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#--------------------------------------------------------------------------------#
#- run settings file
scriptPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
. $scriptPath/scriptSettings.sh $1
#--------------------------------------------------------------------------------#
#- run case settings file
. $casePath/runSettings.sh "DEM"
#--------------------------------------------------------------------------------#

echo  $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

#- call function to run DEM case
parDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

#- return
cd $currentPath
