# adapt particle parameters in deadman region

group DMcoke region DM
set group DMcoke property/atom/polydispparcel effvolfactor 1.5

group nonDMcoke subtract coke DMcoke
set group nonDMcoke property/atom/polydispparcel effvolfactor ${poly_C}

group DMcoke clear
group nonDMcoke clear
