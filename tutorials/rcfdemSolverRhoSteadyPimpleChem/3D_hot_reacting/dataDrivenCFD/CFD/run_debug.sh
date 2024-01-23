mpirun -np 8 valgrind --leak-check=full --show-reachable=yes --track-origins=yes --log-file=nc.vg.%p rcfdemSolverRhoSteadyPimpleChem -parallel
