cd rCFD
cd CFD
rm -r dynamicCode
mpirun -np 4 rcfdemSolverCoupledHeattransfer -parallel > runRCFD.log 2>&1
