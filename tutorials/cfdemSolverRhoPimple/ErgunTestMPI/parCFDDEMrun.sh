#!/bin/bash

cd CFD
blockMesh
decomposePar -force
mpirun -np 4 cfdemSolverRhoPimple -parallel
