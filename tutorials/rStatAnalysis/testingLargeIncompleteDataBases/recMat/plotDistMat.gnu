reset
clear

set terminal cairolatex pdf dash dl 1

set pm3d map
set contour surface
set cntrparam levels discr 10 
set samples 50
set isosamples 50 

set palette maxcolors 8
set palette defined(\
0 0.2314 0.2980 0.7529,\
0.125000 0.384300 0.509800 0.917600,\
0.250000 0.552900 0.690200 0.996100,\
0.375000 0.721600 0.815700 0.976500,\
0.500000 0.866700 0.866700 0.866700,\
0.625000 0.960800 0.768600 0.678400,\
0.750000 0.956900 0.603900 0.482400,\
0.875000 0.870600 0.376500 0.302000,\
1 0.7059 0.0157 0.1490\
)

set cbtics 0.6
set mcbtics 8
#set cblabel 'distance'
#unset cbtics
set cbrange [0:0.6]
set format cb "%.1f"

#set xrange[0:3]
#set yrange[0:3]

set xlabel 't [s]'
set ylabel 't [s]'

# evil manual meddling with axes tics
set xtics ("0" 0, "500" 50, "1000" 100, "1500" 150, "2000" 199)
set ytics ("0" 0, "500" 50, "1000" 100, "1500" 150, "2000" 199)
#set mxtics 2
#set mytics 2

set size square


set out 'myDistMatrix.tex'
splot 'myMatrix.txt' matrix using ($1):($2):(1-$3) with image notitle
