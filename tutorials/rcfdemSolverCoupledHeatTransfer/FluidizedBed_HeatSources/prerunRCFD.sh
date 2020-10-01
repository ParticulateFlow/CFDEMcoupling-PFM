cd rCFD
mv ../db1/DEM/liggghts.restartCFDEM* DEM/liggghts.restart

cd CFD
rm -r dataBase*
rm -r proc*
rm -r 0
cp -r orig.0 0

mv ../../db1/CFD/dataBase1 .
mv ../../db2/CFD/dataBase2 .
decomposePar
./decomposeParDB.sh dataBase1
./decomposeParDB.sh dataBase2
