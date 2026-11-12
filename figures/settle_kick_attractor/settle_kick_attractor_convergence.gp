reset
#set terminal pngcairo enhanced size 1200,850 font "Helvetica,16"
#set output "settle_kick_attractor_convergence_gnuplot.png"

datafile = "settle_kick_attractor_convergence.dat"
dc = 0.5016666666666802

set title "Settle-kick evolution of the kicked-block median delay"
set xlabel "Chorus cycles after kick, n"
set ylabel "Kicked-block median delay, {/Symbol-Oblique d}_{B}(n)"
set xrange [1:21]
set yrange [-3.15:0.95]
set xtics 2
set grid back lc rgb "#dddddd" lw 1
set key bottom right opaque box
set border lw 1.2

set arrow 1 from graph 0, first 0 to graph 1, first 0 nohead lc rgb "#555555" lw 1
set arrow 2 from graph 0, first dc to graph 1, first dc nohead dt 2 lc rgb "#d62728" lw 2
set arrow 3 from graph 0, first -dc to graph 1, first -dc nohead dt 2 lc rgb "#2457d6" lw 2
set label 1 sprintf("+d_c = +%.6f",dc) at 14.7,dc+0.09 tc rgb "#d62728"
set label 2 sprintf("-d_c = -%.6f",dc) at 14.7,-dc-0.12 tc rgb "#2457d6"

plot datafile using 1:3 with linespoints lw 2 pt 7 ps 0.75 lc rgb "#777777" title "zero-kick control", \
     datafile using 1:5 with linespoints lw 2.5 pt 7 ps 0.75 lc rgb "#2457d6" title "negative kick, {/Symbol D}{/Symbol f}=-0.008", \
     datafile using 1:7 with linespoints lw 2.5 pt 7 ps 0.75 lc rgb "#d62728" title "positive kick, {/Symbol D}{/Symbol f}=+0.00575"
