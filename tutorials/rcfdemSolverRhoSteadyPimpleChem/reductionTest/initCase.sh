#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

cd DEM
mpirun -np 8 $CFDEM_LIGGGHTS_BIN_DIR/liggghts -in in.liggghts_init > init.log 2>&1
