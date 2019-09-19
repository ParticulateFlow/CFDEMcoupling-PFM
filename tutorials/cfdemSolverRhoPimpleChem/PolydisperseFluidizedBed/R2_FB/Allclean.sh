#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

cd CFD
cleanCase
cd - 

rm -f log*

rm ./DEM/post/dump*

# ----------------------------------------------------------------- end-of-file
