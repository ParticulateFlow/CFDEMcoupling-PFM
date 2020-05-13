./Allclean.sh
cd system
m4 -P ../../../CFDDEM/CFD/system/blockMeshDict_slot.m4 > blockMeshDict
cd ..
blockMesh
mirrorMesh -overwrite

topoSet -dict system/topoSetDict
createPatch -overwrite
rm -r 0
cp -r orig.0/ 0

setsToZones -noFlipMap

setFields
decomposePar -force
