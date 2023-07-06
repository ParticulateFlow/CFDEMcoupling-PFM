./preCreateDB.sh
cd CFDDEM/CFD
rm -r postProcessing # make sure no old info on cellObj1 exists, otherwise current data is written to different file and cannot be found by updateVelScaling.py

mpirun -np 32 cfdemSolverRhoPimple -parallel > CFDDEM.log 2>&1

cd ../..
