#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest CFD part
# Christoph Goniva - March 2012
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_cfdemSolverPimple_ErgunTestMPI_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPimple" #"cfdemSolverPimpleDyM_22x"#
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

#------------------------------#
#  octave

#- change path
cd octave

#- rmove old graph
rm cfdemSolverPiso_ErgunTestMPI.eps

#- run octave
octave totalPressureDrop.m

#- show plot 
evince cfdemSolverPimple_ErgunTestMPI.eps
#------------------------------#

#- copy log file to test harness
cp ../../$logfileName $testHarnessPath
cp cfdemSolverPimple_ErgunTestMPI.eps $testHarnessPath

#- clean up case
cd ..
rm -rf 0.*
rm -rf processor*
rm -rf patchAverage_pressureDrop
rm -rf probes
rm -rf postProcessing
rm log.liggghts
rm ../DEM/post/dump*
rm -rf particles

#- preserve post directory
echo "dummyfile" >> $casePath/DEM/post/dummy


