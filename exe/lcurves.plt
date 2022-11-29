set terminal png size 800,600 enhanced font "Helvetica,20"
set output 'lightcurves.png'
set style data lines
set title 'Lightcurves'
set xtics font ",18"
set ytics font ",18"
set xlabel "Time [s]" font ",24"
set ylabel "log L [erg/s]" font ",24"
#set xrange [5.9e4:6.4e4]
#set yrange [37.6:39.6]
plot \
'./lum.dat' using 1:2 lt rgb "blue" lw 4 title "Disk", \
'./lum.dat' using 1:3 lt rgb "green" lw 2 title "Corona",\






