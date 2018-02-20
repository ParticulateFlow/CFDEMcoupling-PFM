set terminal qt 0
set title "Fractional Reduction"
set xtic auto
set ytic auto
set xlabel "Time, (s)"
set ylabel "Reduction Rate, (%)"
set grid xtics ytics
set autoscale
plot "Output.dat" using 1:9 title 'fw' with lines lc rgb '#0060ad' lt 2.5 lw 2.5, \
	"Output.dat" using 1:10 title 'fm' with lines lc rgb '#dd181f' lt 2.5 lw 2.5, \
	"Output.dat" using 1:11 title 'fh' with lines lt 2.5 lw 2.5,\
	"Output.dat" using 1:($11*1/9+$10*2/9+$9*6/9) title 'ftot' with lines lc "black" lt 2.5 lw 2.5

pause 1
reread
