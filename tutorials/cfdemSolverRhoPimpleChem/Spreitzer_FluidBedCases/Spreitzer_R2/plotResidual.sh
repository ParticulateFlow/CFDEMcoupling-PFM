#!/bin/bash
# Open and plot the residuals for time/iteration and gas species

casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

cd $casePath
#echo "Load Gnuplot"
#gnuplot > load 'Residuals.plt' &
echo "Residual vs Time"
foamMonitor -l CFD/postProcessing/residuals/0/residuals.dat &
echo "gasResidual vs Time"
foamMonitor -l CFD/postProcessing/gasResidual/0/residuals.dat
wait

