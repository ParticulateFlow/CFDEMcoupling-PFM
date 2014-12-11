#!/bin/bash

# adapt settings for restart run
cp ./CFD/constant/liggghtsCommands_restart ./CFD/constant/liggghtsCommands
cp ./CFD/constant/couplingProperties_restart ./CFD/constant/couplingProperties
cp ./CFD/system/controlDict_restart ./CFD/system/controlDict

