#!/bin/bash
#------------------------------------------------------------------------------
# parCFDDEMrun script for HeatTransferRanzMarshall test case
# run HeatTransferRanzMarshall CFD-DEM
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
headerText="run_parallel_cfdemSolverRhoPimple_HeatTransferRanzMarshall"
logfileName="log_$headerText"
solverName="cfdemSolverRhoPimple"
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
runCleanUp="false"
#------------------------------------------------------------------------------

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

if [ $runCleanUp == "true" ]
    then
        #- clean up case
        echo "deleting data at: $casePath :\n"
        source $WM_PROJECT_DIR/bin/tools/CleanFunctions
        cd $casePath/CFD
        cleanCase
        rm $casePath/DEM/post/*.*
        touch $casePath/DEM/post/.gitignore
fi

echo "done"
