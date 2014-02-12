#!/bin/bash

#===================================================================#
# compile routine for LAMMPS m2m library
# Christoph Goniva - Oct. 2012, DCS Computing GmbH
#===================================================================

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

NOW="$(date +"%Y-%m-%d-%H:%M")"
logDir="log"

cd $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc
mkdir -p $logDir

#--------------------------------------------------------------------------------#
#- define variables
logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
logfileName="log_compileM2Mlib"
headerText="$logfileName""-$NOW"
makeFileName="Makefile.$CFDEM_LIGGGHTS_MAKEFILE_NAME"
libraryPath="$CFDEM_M2MLIB_PATH"
#--------------------------------------------------------------------------------#

compileLMPlib $logpath $logfileName $headerText $makeFileName $libraryPath
