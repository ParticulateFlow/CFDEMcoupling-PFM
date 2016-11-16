#!/bin/bash

grep -n "concentration of CO2" log_Particle_in_Duct > CO2concentration-output
sed -r 's/.{27}//' CO2concentration-output > changeOfCO2
rm CO2concentration-output
sed -r 's/^ = //' changeOfCO2 > changeOfCO2-2
rm changeOfCO2
sed -r 's/^=//' changeOfCO2-2 > changeOfCO2-3
rm changeOfCO2-2


grep -n "concentration O2" log_Particle_in_Duct > O2concentration-file
sed -r 's/.{23}//' O2concentration-file > O2_file_2
rm O2concentration-file
sed -r 's/^ = //' O2_file_2 > O2_output
rm O2_file_2
sed -r 's/^=//' O2_output > O2_output_2
rm O2_output


grep -n "rhogas =" log_Particle_in_Duct > rhogas_file
sed -r 's/.{13}//' rhogas_file > rhogas_file_2
rm rhogas_file
sed -r 's/^ = //' rhogas_file_2 > rhogas_output
rm rhogas_file_2
sed -r 's/^=//' rhogas_output > rhogas_output-2
rm rhogas_output

grep -n "mass of particle = " log_Particle_in_Duct > pmass_file
sed -r 's/.{23}//' pmass_file > pmass_file_2
rm pmass_file
sed -r 's/^ = //' pmass_file_2 > pmass_file_3
rm pmass_file_2
sed -r 's/^=//' pmass_file_3 > pmass_file_4
rm pmass_file_3


grep -n "mass of O2 = " log_Particle_in_Duct > O2mass_file
sed -r 's/.{17}//' O2mass_file > O2mass_file_2
rm O2mass_file
sed -r 's/^ = //' O2mass_file_2 > O2mass_file_3
rm O2mass_file_2
sed -r 's/^=//' O2mass_file_3 > O2mass_file_4
rm O2mass_file_3

grep -n "mass of Co2 = " log_Particle_in_Duct > CO2mass-ouput
sed -r 's/.{18}//' CO2mass-ouput > CO2mass-ouput-2
rm CO2mass-ouput
sed -r 's/^ = //' CO2mass-ouput-2 > CO2mass-ouput-3
rm CO2mass-ouput-2
sed -r 's/^=//' CO2mass-ouput-3 > CO2mass-ouput-4
rm CO2mass-ouput-3