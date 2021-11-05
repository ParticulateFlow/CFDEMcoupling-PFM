#!/bin/bash

#===================================================================#
# post-run script for case
# Tim MJ Nijssen - September 2021
#===================================================================#

#--------------------------------------------------------------------------------#
#- run settings file
scriptPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
. $scriptPath/scriptSettings.sh $1
#--------------------------------------------------------------------------------#

#- make results directories
dateTime=$(date +%y%m%d_%H%M%S)
resultsPath="$casePath"/results_"$runTitle"_"$dateTime"
mkdir $resultsPath
mkdir $resultsPath/DEM
mkdir $resultsPath/DEM/dump_run
mkdir $resultsPath/CFD
mkdir $resultsPath/log
mkdir $resultsPath/VTK

#- reconstruct and convert CFD data
cd $casePath/CFD
reconstructPar -noLagrangian
foamToVTK

#- move CFD files
mv $casePath/CFD/0* $resultsPath/CFD/
mv $casePath/CFD/1* $resultsPath/CFD/
mv $casePath/CFD/2* $resultsPath/CFD/
mv $casePath/CFD/3* $resultsPath/CFD/
mv $casePath/CFD/4* $resultsPath/CFD/
mv $casePath/CFD/5* $resultsPath/CFD/
mv $casePath/CFD/6* $resultsPath/CFD/
mv $casePath/CFD/7* $resultsPath/CFD/
mv $casePath/CFD/8* $resultsPath/CFD/
mv $casePath/CFD/9* $resultsPath/CFD/
mv $casePath/CFD/particleProbes* $resultsPath/CFD/
mv $casePath/CFD/postProcessing* $resultsPath/CFD/
mv $casePath/CFD/dynamicCode* $resultsPath/CFD/

mv $casePath/CFD/VTK/CFD* $resultsPath/VTK/

cp -r $casePath/CFD/org* $resultsPath/CFD/
cp -r $casePath/CFD/constant $resultsPath/CFD/
cp -r $casePath/CFD/system $resultsPath/CFD/
cp -r $casePath/CFD/clockData $resultsPath/CFD/

#- copy DEM data
cp $casePath/DEM/* $resultsPath/DEM/
cp -r $casePath/DEM/mesh $resultsPath/DEM/

#- move and convert DEM init files
if ls $casePath/DEM/post/*init* 1> /dev/null 2>&1; then
    mkdir $resultsPath/DEM/dump_init
    mv $casePath/DEM/post/*init* $resultsPath/DEM/dump_init/

    cd $resultsPath/DEM/dump_init
    python2 $CFDEM_LPP_DIR/lpp.py dump*
    mv $resultsPath/DEM/dump_init/*.vtk $resultsPath/VTK/
fi

#- move and convert DEM run files
mv $casePath/DEM/post/*run* $resultsPath/DEM/dump_run/

cd $resultsPath/DEM/dump_run
python2 $CFDEM_LPP_DIR/lpp.py dump*
mv $resultsPath/DEM/dump_run/*.vtk $resultsPath/VTK/

#- copy log files
cp $casePath/log/* $resultsPath/log/
mv $casePath/slurm* $resultsPath/log/

#- return to path
cd $currentPath
