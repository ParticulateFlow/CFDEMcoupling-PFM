#!/bin/bash

#===================================================================#
# compile routine for CFDEMcoupling + LIGGGHTS, part of CFDEMproject 
# Christoph Goniva - May. 2012, DCS Computing GmbH
#                    update March 2014
#===================================================================#

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

NOW="$(date +"%Y-%m-%d-%H:%M")"
logDir="log"


cd $CFDEM_PROJECT_DIR/etc
mkdir -p $logDir

#================================================================================#
# compile LIGGGHTS src
#================================================================================#
bash $CFDEM_PROJECT_DIR/etc/compileLIGGGHTS.sh

#================================================================================#
# compile LIGGGHTS libraries
#================================================================================#
bash $CFDEM_PROJECT_DIR/etc/compileLIGGGHTS_lib.sh

#================================================================================#
# compile CFDEMcoupling
#================================================================================#
bash $CFDEM_PROJECT_DIR/etc/compileCFDEMcoupling.sh
