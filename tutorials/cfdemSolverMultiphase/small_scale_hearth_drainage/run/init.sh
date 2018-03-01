#!/bin/bash

#===================================================================#
# allrun script for cfdemSolverMultiphase 
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
nrProcs=4;

#fill with the particles
cd $casePath/DEM
mpirun -np $nrProcs liggghts < in.liggghts_fill |& tee $casePath/log_init


#Run 5 seconds of coupled simulation with closed outlet to equilibrate the system
cd $casePath/CFD
cp -r 0.org 0
cp system/controlDict.init system/controlDict
cp constant/couplingProperties.init constant/couplingProperties
setFields
decomposePar -force
mpirun -np $nrProcs cfdemSolverMultiphase -parallel |& tee $casePath/log_equilibrate
reconstructPar -latestTime
mv 5 $casePath/DEM/post/restart/

#switch the outlet boundary conditions to open for the run simulation, easier done manually or by using groovyBC...
perl -0777 -i.original -pe 's/    outlet\n    {\n        type            zeroGradient;/    outlet\n    {\n        type        fixedValue;\n        value       uniform 0;/igs' $casePath/DEM/post/restart/5/p_rgh

perl -0777 -i.original -pe 's/    outlet\n    {\n        type            fixedValue;/    outlet\n    {\n        type        inletOutlet;\n        inletValue        uniform (0 0 0);/igs' $casePath/DEM/post/restart/5/U

