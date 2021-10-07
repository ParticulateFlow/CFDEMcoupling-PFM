#!/bin/bash

#===================================================================#
# clean script for case 
# Tim MJ Nijssen - September 2021
# Based on: Christoph Goniva - Feb. 2011
#===================================================================#

#--------------------------------------------------------------------------------#
#- run settings file
scriptPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
. $scriptPath/scriptSettings.sh $1
#--------------------------------------------------------------------------------#

#- clean up case
rm -r $casePath/*~
rm -r $casePath/#*#
rm -r $casePath/DEM/*~
rm -r $casePath/DEM/#*#
rm -r $casePath/CFD/*~
rm -r $casePath/CFD/#*#
rm -r $casePath/*.e*
rm -r $casePath/*.o*

rm -r $casePath/log*

rm -r $casePath/CFD/0*
rm -r $casePath/CFD/1*
rm -r $casePath/CFD/2*
rm -r $casePath/CFD/3*
rm -r $casePath/CFD/4*
rm -r $casePath/CFD/5*
rm -r $casePath/CFD/6*
rm -r $casePath/CFD/7*
rm -r $casePath/CFD/8*
rm -r $casePath/CFD/9*

rm -r $casePath/CFD/log*

rm -r $casePath/CFD/org.0/*~
rm -r $casePath/CFD/constant/*~
rm -r $casePath/CFD/system/*~
rm -r $casePath/DEM/*~

rm -r $casePath/CFD/processor*
rm -r $casePath/CFD/clockData
rm -r $casePath/CFD/VTK
rm -r $casePath/CFD/particleProbes
rm -r $casePath/CFD/postProcessing

rm -r $casePath/DEM/log*
rm -r $casePath/DEM/post/dump*
rm -r $casePath/DEM/post/*.vtk
rm -r $casePath/DEM/post/*.txt

#- liggghts restart
if [ -f "$casePath/DEM/post/restart/liggghts.restart" ]; then
    echo -n "Clean LIGGGHTS restart file(s)? (y/N)? "
    read -t 10 answer
    if echo "$answer" | grep -iq "^y" ;then
	rm -r $casePath/DEM/post/restart/*
    fi
fi

#- CFD mesh
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo -n "Clean CFD mesh? (y/N)? "
    read -t 10 answer
    if echo "$answer" | grep -iq "^y" ;then
    
	rm -r $casePath/CFD/constant/extendedFeatureEdgeMesh
	rm $casePath/CFD/constant/triSurface/*.eMesh
	rm $casePath/CFD/constant/polyMesh/boundary
	rm $casePath/CFD/constant/polyMesh/faces
	rm $casePath/CFD/constant/polyMesh/neighbour
	rm $casePath/CFD/constant/polyMesh/owner
	rm $casePath/CFD/constant/polyMesh/points
    fi
fi

#- results
if ls $casePath/results* 1> /dev/null 2>&1; then
    echo -n "Clean results? (y/N)? "
    read -t 10 answer
    if echo "$answer" | grep -iq "^y" ;then
	rm -r $casePath/results*
    fi
fi

#- function objects
if [ -d "$casePath/CFD/dynamicCode" ]; then
    echo -n "Clean dynamicCode? (y/N)? "
    read -t 10 answer
    if echo "$answer" | grep -iq "^y" ;then
	rm -r $casePath/CFD/dynamicCode
    fi
fi

cd $currentPath


