# adapt particle parameters in deadman region

group DMcoke region DM
set group DMcoke property/atom/polydispparcel effvolfactor 1.35

group hearthcoke region hearth

group nonDMcoke subtract hearthcoke DMcoke
set group nonDMcoke property/atom/polydispparcel effvolfactor ${poly_C}

group DMcoke clear
group hearthcoke clear
group nonDMcoke clear

# only need to update one reaction fix per group
fix_update chem2
