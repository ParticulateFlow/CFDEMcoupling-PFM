./preCreateDB.sh
cd CFDDEM/CFD

mpirun -np 8 cfdemSolverRhoPimple -parallel > CFDDEM.log 2>&1

cd ../..
