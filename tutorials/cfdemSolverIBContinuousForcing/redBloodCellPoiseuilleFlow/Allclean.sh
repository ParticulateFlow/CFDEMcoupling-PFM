#!/bin/bash
#------------------------------------------------------------------------------
# Allclean script for redBloodCellPoiseuilleFlow test case
# clean up redBloodCellPoiseuilleFlow
# Achuth N. Balachandran Nair - October 2018
#------------------------------------------------------------------------------

#- source CFDEM env vars
. ~/.bashrc

#- source CFDEM functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_cfdemSolverIBContinuousForcing_redBloodCellPoiseuilleFlow"
logfileName="log_$headerText"
solverName="cfdemSolverIBContinuousForcing"
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="false"
postproc="true"
#--------------------------------------------------------------------------------#

#- copy log file to test harness
cp ../../$logfileName $testHarnessPath

#- clean up case
echo "deleting data at: $casePath"
source $WM_PROJECT_DIR/bin/tools/CleanFunctions
cd $casePath/CFD
cleanTimeDirectories
rm -rf processor*
rm -r $casePath/CFD/clockData
rm $casePath/DEM/post/*.*
#rm -r $casePath/DEM/post/restart/*.*
#- preserve post directory
touch $casePath/DEM/post/.gitignore
echo "done"

