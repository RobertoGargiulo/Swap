set terminal png # tikz color # standalone
#set output "figures/figure.png"
set size ratio 0.5
# set title "It would be nice if this were converted into a caption!"
#set xtics 10000 nomirror
#set ytics format "%.1f"
#set xtics (1, 10, 1e2, 1e3, 1e4, 1e5, 1e6) format "10^{%T}"
set margins -1,2,0,0
#set logscale x
#set yrange [0:1.1]
#set xrange [1:1e5]
# set key notitle invert under reverse Left left spacing 2 samplen 0.7
# set arrow 1 filled from graph 0.4, 0.7 to graph 0.6, 0.7
# set label 1 at graph 0.5, 0.75 "$k$" center


filename = 'time_avg_nspin10.steps4000.iterations30.time_step0.5.txt' 

### Plot abs(AVG) or FLUCT at varying 'nspin'
type = 'AVG'
Jint = '0.50'
Vint = '1.00'
hz = '3.00'
iter=30
steps=4000
set yrange [-0.1:1.1]
set xrange [-0.1:6]
set xlabel "V/J"
set ylabel "I"
set output "figures/figure.png"
set title "Imbalance plot at varying parameters. \n L = 10, Steps = 4000, Iterations = 30, Time Step = 0.5"
array low[8]
array high[8]
low[1] = 1; low[2] = 6; low[3] = 11; low[4] = 17; low[5] = 23; low[6] = 28;
low[7] = 34; low[8] = 39
high[1] = 5; high[2] = 10; high[3] = 16; high[4] = 22; high[5] = 27; high[6] = 33;
high[7] = 38; high[8] = 43
do for [i=1:8] {
hz(i)=system("awk '{if(NR==".low[i]."+1) print $3/$1}' time_avg_nspin10.steps4000.iterations30.time_step0.5.txt")
system("echo ".low[i]." ".high[i])
}
plot for [i=1:8] filename every ::low[i]-1::high[i]-1 u ($2/$1):(-$4) w lp title 'hz/J = '.hz(i)



### Plot abs(AVG) vs FLUCT 
nspin = '8'
steps = '4000'
n_iter = '30'
kick = '0.00'
Jint = '0.00'
Vint = '1.00'
hz = '3.01'
#set xrange [0.05:2e2]
#plot  'filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J'.Jint.'_V'.Vint.'_h'.hx.'_hz'.hz.'_no_kick'.kick.'.txt' u 2:(-$1) w l title 'AVG',\
      'filetype.'_FLUCT_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J'.Jint.'_V'.Vint.'_h'.hx.'_hz'.hz.'_no_kick'.kick.'.txt' u 2:1 w l title 'FLUCT'






### Plot single iteration
nspin = '8'
h_coupling = '0.15'
kick = '0.05'
steps = 100
type = 'early'
start = 1
#set yrange [-1.1:1.1]
#unset logscale x
#set xtics 10 format "%0.s"
#set title 'nspin = '.nspin.', h = '.h_coupling.', eps = '.kick.', start = '.sprintf("%0.f",start)
#set output "figures/".type."_steps".steps."_nspin".nspin."_h".h_coupling."_eps".kick.".png"
#plot [1:steps] 'filetype.'_nspin'.nspin.'_steps100000_iterations40_h0.30_kick0.10.txt' every ::start::start+steps u ($2-start):1 w fillsteps title 'iteration = 1',\
#              'filetype.'_nspin'.nspin.'_steps100000_iterations40_h0.30_kick0.10.txt' every ::1e5+1::1e5+steps u 2:1 w fillsteps title 'iteration = 2'
