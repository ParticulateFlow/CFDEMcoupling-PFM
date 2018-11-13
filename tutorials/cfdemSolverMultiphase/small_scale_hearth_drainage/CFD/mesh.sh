#!/bin/bash
m4 constant/polyMesh/blockMeshDict.m4 > constant/polyMesh/blockMeshDict
m4 system/topoSetDict.m4 > system/topoSetDict
blockMesh
topoSet
createPatch -overwrite

