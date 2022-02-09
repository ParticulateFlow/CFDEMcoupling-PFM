#!/bin/bash

#===================================================================#
# allrun script
# Tim MJ Nijssen - September 2021
# based on: Christoph Goniva - August 2011
#===================================================================#

#--------------------------------------------------------------------------------#
#- run settings file
scriptPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
. $scriptPath/scriptSettings.sh $1
#--------------------------------------------------------------------------------#

#- do all tasks before actual run
echo "runAllPar: calling preRunAllPar"
. $scriptPath/preRunAllPar.sh $casePath

# duplicate restart file to restartSequence
echo "runAllPar: duplicating LIGGGHTS restart file"
cp $casePath/DEM/post/restart/liggghts.restart $casePath/DEM/post/restart/liggghts.restartSequence

#- run parallel CFD-DEM
echo "runAllPar: Running CFD-DEM"
. $scriptPath/runCFDDEMPar.sh $casePath

#- do all tasks after actual run
echo "runAllPar: calling postRunAllPar"
. $scriptPath/postRunAllPar.sh $casePath

