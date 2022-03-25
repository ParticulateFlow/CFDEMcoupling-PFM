#!/bin/bash
#------------------------------------------------------------------------------
# allrun script for falling sphere testcase
# run falling_sphere_two_way_coupling
# Achuth N. Balachandran Nair - Oct. 2018
#------------------------------------------------------------------------------

source $CFDEM_PROJECT_DIR/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

echo $casePath

# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd CFD/
    blockMesh
fi

#- run parallel CFD-DEM
bash $casePath/parCFDDEMrun.sh
