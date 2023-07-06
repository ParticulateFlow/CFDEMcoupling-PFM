#./pre_init.sh

#cd initLayerStructure/DEM
#mpirun -np 32 liggghts < in.liggghts_fill #> fill.log 2>&1

cd ../CFD
mpirun -np 32 cfdemSolverRhoPimple -parallel #> createLayers.log 2>&1

cd ../..
./post_init.sh
