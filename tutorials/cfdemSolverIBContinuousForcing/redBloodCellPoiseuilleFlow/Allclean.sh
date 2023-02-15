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

#------------------------------------------------------------------------------
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

#------------------------------------------------------------------------------
#- clean up case
echo "deleting data at: $casePath"
source $WM_PROJECT_DIR/bin/tools/CleanFunctions
cd $casePath/CFD
cleanTimeDirectories
rm -rf processor* > /dev/null 2>&1
rm -rf clockData > /dev/null 2>&1
echo "Cleaning $casePath/DEM/post"
rm $casePath/DEM/post/*.*  > /dev/null 2>&1
#rm -r $casePath/DEM/post/restart/*.*
#touch $casePath/DEM/post/restart/.gitignore
echo "done"

