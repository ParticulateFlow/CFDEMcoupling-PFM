cd CFD
rm -r dynamicCode
rm -r [0-9]*
rm -r proc*
rm -r constant/polyMesh
rm -r postProcessing
m4 system/blockMeshDict.m4 > system/blockMeshDict
blockMesh
topoSet -dict system/topoSetDict

createPatch -overwrite
cp -r orig.0 0

setsToZones -noFlipMap

decomposePar -force

foamToSurface surface.stl
surfaceSplitByPatch surface.stl
rm -r stl_files
mkdir stl_files
mv surface_*.stl stl_files/
