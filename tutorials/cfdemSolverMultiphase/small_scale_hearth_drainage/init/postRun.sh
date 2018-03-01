#!/bin/bash

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

#- get VTK data from liggghts dump file
cd $casePath/DEM/post
python $CFDEM_LPP_DIR/lpp.py dump*.liggghts_run

#- get VTK data from CFD sim
cd $casePath/CFD
reconstructPar -withZero -noLagrangian
foamToVTK

cd $casePath
cp -r ./CFD/1/ ../run/CFD/1.org
cp ./DEM/post/restart/liggghts.restartCFDEM_1.000000 ../run/DEM/post/restart/liggghts.restart

