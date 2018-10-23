#!/bin/bash

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

#- clean up case
echo "deleting data at: $casePath :\n"
source $WM_PROJECT_DIR/bin/tools/CleanFunctions
cd $casePath/CFD
mv 1.org ..
cleanCase
mv ../1.org .
rm -r $casePath/CFD/clockData
rm $casePath/DEM/post/*.*
touch $casePath/DEM/post/.gitignore
rm $casePath/log*
rm $casePath/DEM/log*
#rm -r $casePath/CFD/0

#echo "Remove restart files?"
#echo "Enter: yes, Ctrl + C: no"
#read

#rm $casePath/DEM/post/restart/*.*


