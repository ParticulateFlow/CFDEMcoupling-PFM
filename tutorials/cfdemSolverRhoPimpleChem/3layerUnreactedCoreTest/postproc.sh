#!/bin/bash

grep -n "x0_:" log_Particle_in_Duct > x0_values
grep -n "check N" log_Particle_in_Duct > total_mole
grep -n "check mass frac" log_Particle_in_Duct > mass_frac
grep -n "dY_" log_Particle_in_Duct > delta_reduction_rate
grep -n "dmA_" log_Particle_in_Duct > layer_mass_transfer
grep -n "dens_" log_Particle_in_Duct > layer_Densities
grep -n "pre-particle density" log_Particle_in_Duct > total_particle_density
grep -n "pre-particle mass" log_Particle_in_Duct > pre_particle_mass
grep -n "post-particle mass" log_Particle_in_Duct > post_particle_mass
grep -n "post redox radius of" log_Particle_in_Duct > post_particle_radius


