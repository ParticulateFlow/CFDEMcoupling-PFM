# adapt particle parameters in deadman region

fix eeDM1 all execute ${Nset} " group DMcoke region DM "
fix eeDM2 all execute ${Nset} " set group DMcoke property/atom/polydispparcel effvolfactor 1.5 "

fix eeDM3 all execute ${Nset} " group nonDMcoke subtract all DMcoke "
fix eeDM4 all execute ${Nset} " set group nonDMcoke property/atom/polydispparcel effvolfactor ${poly_C} "

fix eeDM5 all execute ${Nset} " group DMcoke clear "
fix eeDM6 all execute ${Nset} " group nonDMcoke clear "
