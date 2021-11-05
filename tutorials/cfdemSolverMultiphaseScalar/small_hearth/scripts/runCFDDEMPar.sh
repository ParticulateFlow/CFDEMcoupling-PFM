#!/bin/bash

#===================================================================#
# CFDDEMrun script for case
# Tim MJ Nijssen - September 2021
# based on: Christoph Goniva - May. 2011
#===================================================================#

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#--------------------------------------------------------------------------------#
#- run settings file
scriptPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
. $scriptPath/scriptSettings.sh $1
#--------------------------------------------------------------------------------#
#- run case settings file
. $casePath/runSettings.sh "CFDDEM"
#--------------------------------------------------------------------------------#

cd $casePath/CFD

    if [ $debugMode == "on" ]; then
        debugMode="valgrind"
    elif [ $debugMode == "strict" ]; then
        debugMode="valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes"  
    else
        debugMode=""
    fi

#- make proc dirs visible
count=0
for i in `seq $nrProcs`
do
    let count=$i-1
    (cd $casePath/CFD/processor$count && touch file.foam)
done

#- header
echo 2>&1 | tee -a /$logpath/$logfileName
echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
echo 2>&1 | tee -a $logpath/$logfileName

#- write path
pwd 2>&1 | tee -a $logpath/$logfileName
echo 2>&1 | tee -a $logpath/$logfileName

#- run applictaion
if [[ $machineFileName == "none" ]]; then
    mpirun -np $nrProcs $debugMode $solverName -parallel 2>&1 | tee -a $logpath/$logfileName
else
    mpirun -machinefile $machineFileName -np $nrProcs $debugMode $solverName -parallel 2>&1 | tee -a $logpath/$logfileName
fi

#- return
cd $currentPath
