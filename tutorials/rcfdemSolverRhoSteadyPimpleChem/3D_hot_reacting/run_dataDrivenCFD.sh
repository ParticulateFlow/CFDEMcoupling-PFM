./preDataDrivenCFD.sh
cd dataDrivenCFD/CFD

mpirun -np 8 topoSet -dict system/topoSetDictCZ -parallel > topoSet.log 2>&1
mpirun -np 8 rcfdemSolverRhoSteadyPimpleChem -parallel > dataDrivenCFD.log 2>&1

cd ../..
