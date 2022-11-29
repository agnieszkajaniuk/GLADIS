set terminal png size 800,600 enhanced font "Helvetica,20"
set output 'sigte_time.png'
set style data lines
set title 'Stability tracks'
set logscale x
set logscale y
set xtics font ",18"
set ytics font ",18"
set xlabel "log Sigma [g/s]" font ",24"
set ylabel "log Teff [K]" font ",24"
set xrange [1.e3:8.e6]
set yrange [5.e5:8.e7]
plot \
'./sigte.020' using 1:2  with points lt rgb "blue" lw 4 title "time evol", \
'./sigma.020' using 1:2 lt rgb "magenta" lw 4 title "Stability curve, R=10.0", \




