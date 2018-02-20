#!/bin/bash

#grep "dmB\[0]" log_test_case > dmB0.dat
#grep "dmB\[1]" log_test_case > dmB1.dat
#grep "dmB\[2]" log_test_case > dmB2.dat
#grep "dmB\[3]" log_test_case > dmB3.dat
#grep -n "pre-layerMass[0]:" log_3layerUnreactedShrinkingCore > pre_particle_mass0.dat
#grep -n "pre-layerMass[1]:" log_3layerUnreactedShrinkingCore > pre_particle_mass1.dat
#grep -n "pre-layerMass[2]:" log_3layerUnreactedShrinkingCore > pre_particle_mass2.dat
#grep -n "pre-layerMass[3]:" log_3layerUnreactedShrinkingCore > pre_particle_mass3.dat
grep -n "post-layerMass[0]:" log_test_case > post_particle_mass0.dat
grep -n "post-layerMass[1]:" log_test_case > post_particle_mass1.dat
grep -n "post-layerMass[2]:" log_test_case > post_particle_mass2.dat
grep -n "post-layerMass[3]:" log_test_case > post_particle_mass3.dat
#grep -n "x0_eq :" log_3layerUnreactedShrinkingCore > x0_eq_values.dat
#grep -n "x0_:" log_3layerUnreactedShrinkingCore > x0_values.dat
#grep -n "dY_" log_3layerUnreactedShrinkingCore > delta_reduction_rate.dat
#grep -n "dmA_" log_3layerUnreactedShrinkingCore > layer_mass_transfer.dat
