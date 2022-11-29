set terminal png size 800,600 enhanced font "Helvetica,20"
set output 'Sigma_time.png'
set style data lines
set title 'Disk surface density'
set logscale x
set logscale y
set xtics font ",18"
set ytics font ",18"
set xlabel "log radius [r_g]" font ",24"
set ylabel "log Sigma [g/s]" font ",24"
set xrange [1:1.e3]
set yrange [1.e3:8.e6]
plot \
'./rsig.2.000000E+04' using 1:2 lt rgb "blue" lw 4 title "t=2.0e4 s", \
'./rsig.2.100000E+04' using 1:2 lt rgb "green" lw 2 title "t=2.1e4 s",\
'./rsig.2.200000E+04' using 1:2 lt rgb "magenta" lw 4 title "t=2.2e4 s", \
'./rsig.2.300000E+04' using 1:2 lt rgb "red" lw 2 title "t=2.3e4 s",\



