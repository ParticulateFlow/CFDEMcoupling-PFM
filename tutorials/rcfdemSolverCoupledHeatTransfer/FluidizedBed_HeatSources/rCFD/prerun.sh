mv ../db1/DEM/liggghts.restartCFDEM_-0.005000 DEM/liggghts.restart

cd CFD

mv ../../db1/CFD/dataBase1 .
mv ../../db2/CFD/dataBase2 .
decomposePar
./decomposeParDB.sh dataBase1
./decomposeParDB.sh dataBase2
