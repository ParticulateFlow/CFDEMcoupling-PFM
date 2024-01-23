cd initLayerStructure/CFD
rm -r 0

cp -r ../../CFDDEM/CFD/orig.0 0
cp 0/partTemp 0/partTempRef

./makePointsInCZ.sh
setFields
decomposePar -force

cd ../..
