cd CFD
rm -r dynamicCode
mpirun -np 4 rcfdemSolverCoupledHeattransfer -parallel | tee run.log
