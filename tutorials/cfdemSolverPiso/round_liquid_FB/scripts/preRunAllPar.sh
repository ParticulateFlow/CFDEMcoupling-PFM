#!/bin/bash

#===================================================================#
# all tasks before actual CFD-DEM run
# Tim MJ Nijssen - September 2021
# based on: Christoph Goniva - August 2011
#===================================================================#

#--------------------------------------------------------------------------------#
#- run settings file
scriptPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
. $scriptPath/scriptSettings.sh $1
#--------------------------------------------------------------------------------#

runTitle=$(basename "$casePath")
#- prompt for costum title
echo -n "Add custom simulation title? (y/N)?"
read -t 10 answer
if echo "$answer" | grep -iq "^y" ;then
    echo -n "Enter simulation title:"
    read runTitle
fi

#- run DEM init
if [ -f "$casePath/DEM/in.liggghts_init" ]; then
    if [ -f "$casePath/DEM/post/restart/liggghts.restart" ];  then
	echo "preRunAllPar: Using existing restart file"
    else
	echo "preRunAllPar: Running DEM init"
	. $scriptPath/runDEMPar.sh $casePath
    fi
fi

#- copy 0.org
echo "preRunAllPar: Copying 0.org"
cp -r $casePath/CFD/org.0 $casePath/CFD/0

#- mesh
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "preRunAllPar: using old mesh"
else
    echo "preRunAllPar: Building mesh"
    . $scriptPath/buildMesh.sh $casePath > $logpath/log.mesh
fi

#-  set fields
if [ -f "$casePath/CFD/system/setFieldsDict" ]; then
    echo "preRunAllPar: Setting fields"
    cd $casePath/CFD
    setFields
fi

#-  decompose
echo "preRunAllPar: Decomposing"
cd $casePath/CFD
decomposePar
cd $currentPath
