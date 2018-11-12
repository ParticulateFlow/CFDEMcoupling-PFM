cd CFD
cp system/controlDict_equil system/controlDict
decomposePar -force
mpirun -np 4 cfdemSolverRhoPimple -parallel
cp system/controlDict_record system/controlDict
mpirun -np 4 cfdemSolverRhoPimple -parallel
reconstructPar

# make sure that time 0 folder is called 0
mv ./*e-1[0-9]* 0

# remove uniform folder
find . -name \uniform -type d -exec rm -rf {} \;

# make database
mkdir dataBase1
mv [0-2]* dataBase1/

rm -rf proc*
rm -rf ./-0.*
