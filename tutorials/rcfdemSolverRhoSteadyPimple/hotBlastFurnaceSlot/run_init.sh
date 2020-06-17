cd initLayerStructure/DEM

mpirun -np 8 liggghts < in.liggghts_fill > fill.log 2>&1
cd ../CFD
cp -r orig.0 0
setFields
decomposePar
mpirun -np 8 cfdemSolverRhoPimple -parallel > createLayers.log 2>&1

cd ../..
./post_init.sh
