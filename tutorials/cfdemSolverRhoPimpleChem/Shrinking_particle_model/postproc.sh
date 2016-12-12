#!/bin/bash

grep -n "Co2 Mass" log_Particle_in_Duct > CO2concentration-output
sed -r 's/.{12}//' CO2concentration-output > changeOfCO2
rm CO2concentration-output
sed -r 's/^s //' changeOfCO2 > changeOfCO2-2
rm changeOfCO2
sed -r 's/^ss //' changeOfCO2-2 > changeOfCO2-3
rm changeOfCO2-2
sed -r 's/^ass //' changeOfCO2-3 > changeOfCO2
rm changeOfCO2-3


grep -n "O2 Mass" log_Particle_in_Duct > O2concentration-file
sed -r 's/.{12}//' O2concentration-file > O2_file_2
rm O2concentration-file
sed -r 's/^s //' O2_file_2 > O2_output
rm O2_file_2
sed -r 's/^ss //' O2_output > O2_output_2
rm O2_output
sed -r 's/^ass //' O2_output_2 > changeOfO2
rm O2_output_2


grep -n "Gas Density" log_Particle_in_Duct > rhogas_file
sed -r 's/.{16}//' rhogas_file > rhogas_file_2
rm rhogas_file
sed -r 's/^y //' rhogas_file_2 > rhogas_output
rm rhogas_file_2
sed -r 's/^ty //' rhogas_output > rhogas_
rm rhogas_output

grep -n "Particle Mass" log_Particle_in_Duct > pmass_file
sed -r 's/.{18}//' pmass_file > pmass_file_2
rm pmass_file
sed -r 's/^s //' pmass_file_2 > pmass_file_3
rm pmass_file_2
sed -r 's/^ss //' pmass_file_3 > pmass_file_4
rm pmass_file_3
sed -r 's/^ass //' pmass_file_4 > pmass_
rm pmass_file_4


