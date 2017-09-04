#!/bin/bash

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
export casePath

echo "generating Residuals"
$casePath
gnuplot Residuals -