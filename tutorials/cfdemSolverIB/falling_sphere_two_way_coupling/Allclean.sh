#!/bin/bash
#------------------------------------------------------------------------------
# Allclean script for falling sphere test case
# clean CFD and DEM part
# Achuth N. Balachandran Nair - October 2018
#------------------------------------------------------------------------------

#- source environment variables
. ~/.bashrc

#- source CFDEM functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

#- clean up case
echo "deleting data at: $casePath"
source $WM_PROJECT_DIR/bin/tools/CleanFunctions
cd $casePath/CFD
cleanCase
rm -r $casePath/CFD/clockData
rm $casePath/DEM/post/*.*
#- preserve post directory
touch $casePath/DEM/post/.gitignore
echo "done"

