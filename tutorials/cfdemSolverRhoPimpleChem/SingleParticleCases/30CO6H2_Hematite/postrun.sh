#!/bin/bash
#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

#- get VTK data from liggghts dump file
cd $casePath/DEM/post
python -i $CFDEM_LPP_DIR/lpp.py dump*.liggghts_run

#- get VTK data from CFD sim
cd $casePath/CFD
reconstructPar
foamToVTK


