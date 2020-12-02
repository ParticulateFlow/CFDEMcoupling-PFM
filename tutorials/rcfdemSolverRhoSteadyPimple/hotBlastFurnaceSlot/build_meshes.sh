# initLayerStructure
cd initLayerStructure/CFD
bash mesh_slot.sh > mesh.log 2>&1
cd ../..

# CFDDEM
cd CFDDEM/CFD
bash mesh_slot.sh > mesh.log 2>&1
cd ../..

# dataDrivenCFD
cd dataDrivenCFD/CFD
bash mesh_slot.sh > mesh.log 2>&1
cd ../..
