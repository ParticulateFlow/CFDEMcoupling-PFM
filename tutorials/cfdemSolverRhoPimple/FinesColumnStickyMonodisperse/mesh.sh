cd CFD
rm -r 0
rm -r constant/polyMesh
blockMesh
topoSet
createPatch -overwrite
setsToZones -noFlipMap
cp -r orig.0 0
decomposePar
