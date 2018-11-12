mv ../db1/DEM/liggghts.restartCFDEM_-0.005000 DEM/liggghts.restart

cd CFD
cp ../../db1/CFD/0/phiMean 0/phiRec
cp ../../db1/CFD/0/pMean 0/pRec
cp ../../db1/CFD/0/rhoMean 0/rhoRec
cp ../../db1/CFD/0/UMean 0/URec
cp ../../db1/CFD/0/UsMean 0/UsRec
cp ../../db1/CFD/0/voidfractionMean 0/voidfractionRec

mv ../../db1/CFD/dataBase1 .
mv ../../db2/CFD/dataBase2 .
decomposePar
./decomposeParDB.sh dataBase1
./decomposeParDB.sh dataBase2
