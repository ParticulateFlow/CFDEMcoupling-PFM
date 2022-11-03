#===================================================================#
# small_hearth
# Tim MJ Nijssen - Oktober 2021
#===================================================================#

This tutorial case provides a demonstration of the extended liquid-solid coupling by Nijssen et al. (2020). The case represents a liquid-solid fluidised bed pellet softening reactor, as described by Nijssen et al. (2021). 

#===================================================================#

The following scripts are available in the /scripts directory:

- cleanCase.sh: Cleans the case directory from all previous output. Prompts the user whether to remove LIGGGHTS restart files, CFD mesh files and result directories
- runAllPar.sh: Runs LIGGGHTS initialisation if no previous restart file is found, then performs a coupled run according to the start and end times set in CFD/system.controlDict. Afterwards, the case is reconstructed and results moved to a separate directory.

#===================================================================#

T.M.J. Nijssen, J.A.M. Kuipers, J. van der Stel, A.T. Adema, K.A. Buist. Complete liquid-solid momentum coupling for unresolved CFD-DEM simulations. International Journal of Multiphase Flow, 2020.
T.M.J. Nijssen, O.J.I. Kramer, P.J. de Moel, J. Rahman, J.P. Kroon, P. Berhanu, E.S. Boek, K.A. Buist, J.P. van der Hoek, J.T. Padding, J.A.M. Kuipers. Experimental and numerical insights into heterogeneous liquid-solid behaviour in drinking water softening reactors. Chemical Engineering Science: X, 2021.
