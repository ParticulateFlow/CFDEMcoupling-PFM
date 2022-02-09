#!/bin/bash

#===================================================================#
# post run script for testcase
# postprocess falling_sphere_two_way_coupling
# Achuth N. Balachandran Nair - Oct. 2018
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
#--------------------------------------------------------------------------------#

#- get VTK data from liggghts dump file
cd $casePath/DEM/post
lpp dump*

#- get VTK data from CFD sim
cd $casePath/CFD
reconstructParMesh -constant -mergeTol 1e-06
reconstructPar -noLagrangian
foamToVTK

#- start paraview
paraview

