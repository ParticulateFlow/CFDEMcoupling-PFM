#!/bin/bash

#===================================================================#
# sets variables for CFDEM utility scripts
# Tim MJ Nijssen - September 2021
#===================================================================#

#- casepath

currentPath=$(pwd)

if [ -z "$1" ]; then
    #- parent directory if no other path is specified
    casePath="$(dirname "$scriptPath")"
else
    #- read input
    casePath="$(dirname "$(readlink -f $1)/.")"
fi

#- other paths
logPath=$casePath/log
mkdir $logPath

