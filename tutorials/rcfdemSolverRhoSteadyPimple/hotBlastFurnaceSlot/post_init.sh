# this script assumes that the initial layer structure has already been created in initLayerStructure
# furthermore, the meshes need to exist

echo "0" > runIndex
cp initLayerStructure/CFD/0/partTemp dataDrivenCFD/CFD/orig.0/partTemp
cp initLayerStructure/CFD/0/partTempRef dataDrivenCFD/CFD/orig.0/partTempRef
cp initLayerStructure/DEM/liggghts.restartCFDEM* CFDDEM/DEM/liggghts.restart
cd CFDDEM/CFD
cp -r orig.0 0
cd ../../dataDrivenCFD/CFD
cp -r orig.0 0
cd ../DEM
echo "variable latestTimeStep equal 0" > latestTimeStep
