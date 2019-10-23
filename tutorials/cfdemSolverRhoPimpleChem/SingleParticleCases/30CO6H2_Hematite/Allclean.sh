#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

cd CFD
cleanCase
cd - 

rm -f log*

rm -r ./DEM/post
mkdir ./DEM/post
mkdir ./DEM/post/restart
touch ./DEM/post/restart/.gitignore

# ----------------------------------------------------------------- end-of-file
