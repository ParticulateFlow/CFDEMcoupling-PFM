#===================================================================#
# small_hearth
# Tim MJ Nijssen - September 2021
#===================================================================#

This tutorial case provides a demonstration of various functions:

- The cfdemSolverMultiphaseScalar solver
- The ParmarBassetForce model
- The wallHeatTransferYagi model
- The massTransferGunn model
- The couple/cfd/dissolve command
- using scripted boundary conditions to obtain frictional outflow conditions
- using external scripting to interface with simulations and adapt conditions and settings

#===================================================================#

The following scripts are available in the /scripts directory:

- cleanCase.sh: Cleans the case directory from all previous output. Prompts the user whether to remove LIGGGHTS restart files, CFD mesh files and result directories
- runAllPar.sh: Runs LIGGGHTS initialisation if no previous restart file is found, then performs a coupled run according to the start and end times set in CFD/system.controlDict. Afterwards, the case is reconstructed and results moved to a separate directory.
- runAllParSequence.sh: Runs LIGGGHTS initialisation if no previous restart file is found, then performs a coupled run according to the settings in the sequenceSettings file. This allows for time-based switching of boundary conditions. Afterwards, the case is reconstructed and results moved to a separate directory.
- runAllParTapping.sh: Runs LIGGGHTS initialisation if no previous restart file is found, then performs a coupled run according to the settings in the tappingSettings file. This allows for conditional switching of boundary conditions. Afterwards, the case is reconstructed and results moved to a separate directory.

NOTE: the latter to scripts rely on Perl (www.perl.org) to function. Perl comes included in most linux distibutions, or can be installed through the package manager.
