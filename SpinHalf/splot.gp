set pm3d map
splot 'datafile.dat'
set terminal png crop
set output 'figures/figure.png'
#set cbrange [0:1]
#set xrange [0:3]
#set yrange [0:3]
unset key
replot
unset output
