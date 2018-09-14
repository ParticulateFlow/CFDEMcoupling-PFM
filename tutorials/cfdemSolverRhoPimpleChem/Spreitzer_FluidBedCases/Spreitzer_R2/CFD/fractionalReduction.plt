set terminal qt 0
set title "Fractional Reduction No CG"
set xtic auto
set ytic auto
set xlabel "Time, (s)"
set ylabel "Reduction Rate, (%)"
set grid xtics ytics
set autoscale
plot "Output.dat" using 1:($2/28545)*100 title 'fw' with lines lc rgb '#0060ad' lt 2.5 lw 2.5, \
	"Output.dat" using 1:($3/28545)*100 title 'fm' with lines lc rgb '#dd181f' lt 2.5 lw 2.5, \
	"Output.dat" using 1:($4/28545)*100 title 'fh' with lines lt 2.5 lw 2.5,\
	"Output.dat" using 1:(($4/28545)*1/9+($3/28545)*2/9+($2/28545)*6/9)*100 title 'ftot' with lines lc "black" lt 2.5 lw 2.5

pause 1
reread
