#!/bin/bash

pushd CFD

transformPoints -scale "(0.001 0.001 0.001)"
checkMesh
#- make the linear system more diagonal dominant to speed-up the linear solvers
renumberMesh -overwrite -noFunctionObjects

popd

