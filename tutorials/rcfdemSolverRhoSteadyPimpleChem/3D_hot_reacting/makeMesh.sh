cd CFDDEM/CFD
./Allclean.sh
rm -r constant/polyMesh
rm -r 0
cp -r orig.0 0
python genBlockMeshExtraHigh.py
blockMesh

# make surface mesh for DEM before refinement
foamToSurface surface.stl
surfaceSplitByPatch surface.stl
mkdir ../DEM/stl_files/
mv surface_*.stl ../DEM/stl_files/

rm -r constant/polyMesh
python genBlockMesh.py
blockMesh

# generate RW regions for refinement
python genRWBlockRegs.py 0.3
topoSet -dict system/topoSetDict_RWReg

refineHexMesh RWRegs -overwrite
python genRWBlockRegs.py 0.16
topoSet -dict system/topoSetDict_RWReg

refineHexMesh RWRegs -overwrite
python genRWBlockRegs.py 0.08
topoSet -dict system/topoSetDict_RWReg

topoSet -dict system/topoSetDict_OtherRegs

topoSet -dict system/topoSetDict_RWFaces
createPatch -overwrite

setsToZones -noFlipMap

rm -r 0
cp -r orig.0 0

cd ../..
rm -r dataDrivenCFD/CFD/constant/polyMesh
cp -r CFDDEM/CFD/constant/polyMesh dataDrivenCFD/CFD/constant/
rm -r initLayerStructure/CFD/constant/polyMesh
cp -r CFDDEM/CFD/constant/polyMesh initLayerStructure/CFD/constant/
