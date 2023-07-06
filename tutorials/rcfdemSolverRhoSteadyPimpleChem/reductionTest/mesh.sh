cd CFD/
rm -r 0
rm -r constant/polyMesh
m4 system/blockMeshDict.m4 > system/blockMeshDict
blockMesh
cp -r orig.0 0
decomposePar
./decomposeParDB.sh
