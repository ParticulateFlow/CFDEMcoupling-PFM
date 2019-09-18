#!/bin/sh
#cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
cd $casePath/CFD
cleanCase
rm $casePath/log*

rm $casePath/DEM/post/liggghts_run*
rm $casePath/DEM/post/dump.liggghts_run
#mkdir $casePath/DEM/post
#mkdir $casePath/DEM/post/restart
#cd $casePath/DEM/post/restart 
#touch liggghts.restart

# ----------------------------------------------------------------- end-of-file
