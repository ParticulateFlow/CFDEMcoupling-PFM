#!/bin/bash

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

#- clean up case
echo "deleting data at: $casePath :\n"
source $WM_PROJECT_DIR/bin/tools/CleanFunctions
cd $casePath/CFD
cleanCase
rm -r $casePath/CFD/clockData
rm $casePath/DEM/post/*.*
touch $casePath/DEM/post/.gitignore
rm -r $casePath/CFD/0
rm $casePath/log*
rm $casePath/*.png

echo "Remove restart file?"
echo "Enter: yes, Ctrl + C: no"
read
rm $casePath/DEM/post/restart/*.*
rm $casePath/DEM/post/restart/liggghts.restartCFDEM*

