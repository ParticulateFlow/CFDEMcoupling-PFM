set terminal qt 1
set title "Fractional Reduction"
set xtic auto
set ytic auto
set xlabel "Time, (s)"
set ylabel "Reduction Rate, (%)"
set grid xtics ytics
set autoscale
plot "MassLayers.dat" using 1:2 title 'mL_Fe' with lines lc rgb '#0060ad' lt 2.5 lw 2.5, \
	"MassLayers.dat" using 1:3 title 'mL_w' with lines lc rgb '#dd181f' lt 2.5 lw 2.5, \
	"MassLayers.dat" using 1:4 title 'mL_m' with lines lt 2.5 lw 2.5, \
	"MassLayers.dat" using 1:5 title 'mL_h' with lines lc "black" lt 2.5 lw 2.5

pause 1
reread
