#!/bin/bash
#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

cd $casePath
grep -E 'pdensity' log_ISO4695 > pdensity.txt
grep -E 'porosity' log_ISO4695 > porosity.txt
grep -E 'active layers' log_ISO4695 > activeLayers.txt
grep -E 'rhoeff_' log_ISO4695 > rhoeff.txt
grep -E 'pmass' log_ISO4695 > pmass.txt
grep -E 'pdensity after mass reduction =' log_ISO4695 > densityARed.txt
