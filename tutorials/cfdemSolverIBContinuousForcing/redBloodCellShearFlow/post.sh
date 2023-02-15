#!/bin/bash
#------------------------------------------------------------------------------
# post script for redBloodCellShearFlow test case
# run post-processing
# Achuth N. Balachandran Nair - October 2018
#------------------------------------------------------------------------------

#- source CFDEM env vars
source ~/.bashrc

#- source CFDEM functions
shopt -s expand_aliases
source $CFDEM_PROJECT_DIR/etc/bashrc

#------------------------------------------------------------------------------
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

#------------------------------------------------------------------------------
#- get VTK data from liggghts dump file
cd $casePath/DEM/post
lpp dump*

#- get VTK data from CFD sim
cd $casePath/CFD
reconstructParMesh -constant -mergeTol 1e-06
reconstructPar -noLagrangian

#- start paraview
paraview

