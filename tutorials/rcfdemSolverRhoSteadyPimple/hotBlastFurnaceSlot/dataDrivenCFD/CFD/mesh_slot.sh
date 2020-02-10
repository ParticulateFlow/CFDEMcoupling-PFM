./Allclean
cd constant/polyMesh
m4 -P blockMeshDict_slot.m4 > blockMeshDict
cd ../..
blockMesh
mirrorMesh

topoSet -dict system/topoSetDict
createPatch -overwrite
rm -r 0
cp -r orig.0/ 0

#foamToSurface surface.stl
#surfaceSplitByPatch surface.stl
#mv surface_*.stl stl_files/
