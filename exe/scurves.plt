set terminal png size 800,600 enhanced font "Helvetica,20"
set output 'scurves.png'
set style data lines
set title 'Stability curves'
set logscale x
set logscale y
set xtics font ",18"
set ytics font ",18"
set xlabel "log Sigma [g/s]" font ",24"
set ylabel "log Teff [K]" font ",24"
set xrange [1.e3:8.e6]
set yrange [5.e4:8.e7]
plot \
'./sigma.020' using 1:2 lt rgb "blue" lw 4 title "point 19, R=10.0", \
'./sigma.059' using 1:2 lt rgb "green" lw 2 title "point 59, R=35.6",\
'./sigma.115' using 1:2 lt rgb "magenta" lw 4 title "point 99, R=99.9", \
'./sigma.199' using 1:2 lt rgb "red" lw 2 title "point 199, R=256.8",\



