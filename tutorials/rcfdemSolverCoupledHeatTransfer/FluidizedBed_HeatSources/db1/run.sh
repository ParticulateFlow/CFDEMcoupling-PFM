cd CFD
cp system/controlDict_equil system/controlDict
decomposePar -force
mpirun -np 4 cfdemSolverRhoPimple -parallel
cp system/controlDict_record system/controlDict
mpirun -np 4 cfdemSolverRhoPimple -parallel
