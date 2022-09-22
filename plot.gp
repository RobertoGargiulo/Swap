set terminal png # tikz color # standalone
set size ratio 0.5
# set title "It would be nice if this were converted into a caption!"
set xtics 10000 nomirror
set ytics format "%.1f"
set xtics (1, 10, 1e2, 1e3, 1e4, 1e5, 1e6) format "10^{%T}"
set margins -1,2,0,0
set logscale x
set yrange [0:1.1]
set xrange [1:1e5]
# set key notitle invert under reverse Left left spacing 2 samplen 0.7
# set arrow 1 filled from graph 0.4, 0.7 to graph 0.6, 0.7


# filetype = 'Swap_Sz'
filetype = 'Clean_MBL_Imbalance'
filetype = 'data/magnetizations/Clean_MBL_OMP_AVG_FLUCT_Imbalance_'

### Plot abs(AVG) or FLUCT at varying 'nspin'
filetype="data/magnetizations/Sz0_DENSE_SWAP_hz_Disorder_AVG_FLUCT_Imbalance_nspin"

nspin = "6"
J = "0.01"
V = "0.25"
hz = "6.00"
steps = "4000"
period = "1.00"
iter = "100"
kick="0.00"

density=1

set logscale x
set yrange [-0.1:1.1]
set xrange [period:period*steps]
#unset xtics
set ylabel "Imbalance"
set xlabel "Time"

nspin = "8"
set title "n_{iter} = ".iter.", period = ".period.", nspin = ".nspin.", V = ".V.", hz = ".hz
#set output "figures/Avg_Imbalance_Evolution_Sz0_SWAP_hz_Disorder_period".period."_V".V."_hz".hz."_kick".kick."_up_to_J0.1.png"
#plot for [J in "0.00 0.02 0.04 0.06 0.08 0.10"] filetype."".nspin."_steps".steps."_period".period."_iterations".iter."_J".J."_V".V."_hz".hz.".txt" every density u 1:(abs($2)) w l title "Swap, J =".J 

param = "nspin16_steps100_period1.00_iterations1_J0.05_V0.50_hz6.00_kick0.05"
file = "data/magnetizations/Sz0_DENSE_SWAP_hz_V_Disorder_Neel_AVG_FLUCT_Imbalance_".param.".txt"
#set terminal tikz color standalone font ",15"
#set output "figures/Single_Run_Imbalance_Evolution_Sz0_SWAP_hz_V_Disorder_".param.".tex"
unset logscale x
set yrange [-1.1:1.1]
set xrange [1:20]
set xtics 5, 5, 50 format "%.0f"
set arrow from 40,0 to 40,100 nohead
y0 = -1.0
y1 = 1.0
do for [x = 1:50] {
    set arrow from x,y0 to x,y1 nohead
}
#plot file u 1:2 w l notitle

kick="0.00"
J = "0.10"
V = "0.25"
list_L = "2 4 6 8 10 12 14 16"
list_iter1 = "100 100 100 100 100 100 100 100"
list_iter2 = "1280 640 320 160 80 80 20 10"
list_iter3 = "2560 1280 640 320 160 80 40 20 10"

list_J1 = "0.00 0.02 0.04 0.06 0.08 0.10"
list_V1 = "0.00 0.25 0.50"
list_hz1 = "0.00 3.00 6.00"

list_J2 = "0.00 0.10 0.20 0.30 0.40 0.50"
list_V2 = "0.00 0.50 1.00 1.50 2.00 2.50 3.00"
list_hz2 = "0.00 1.00 2.00 3.00 4.00 5.00 6.00"
list_kick = "0.00 0.05 0.10 0.15 0.20"

start=0.00; end=4.00; inc=0.10
list(start,end,increment)=system(sprintf("seq %.2f %.2f %.2f", start, increment, end))
list_V3 = list(0.00,4.00,0.10)  #at J = 0.10, hz = 6.00, kick = 0.00
#list_V3 = ""; a=start-inc; while (a<end) {list_V3=list_V3.sprintf(" %.2f",a=a+inc)}
#print list_V3


#Varying L

period = "1.00"
J = "0.05"
V = "0.50"
hz = "6.00"
kick = "0.00"
steps = "10000"

list_L = "2 4 6 8 10 12 14 16"
list_iter = "2560 1280 640 320 160 80 40 20 10"

init_state = "Nayak"
filetype = "data/magnetizations/Sz0_DENSE_SWAP_hz_V_Disorder_".init_state."_AVG_FLUCT_Imbalance_nspin"
figure = "figures/Avg_Imbalance_Evolution_Sz0_SWAP_hz_V_Disorder_".init_state."_Init_State_period".period."_J".J."_V".V."_hz".hz."_kick".kick.".png"

set key left bottom box 3
set xrange [1:steps]
set title "n_{iter} = 2^{1-L/2} 2560, period = ".period.", J = ".J.", V = ".V.", hz = ".hz.", kick = ".kick."\n{/*0.8 ".init_state." initial state}"
#set output figure
#plot for [i=1:6] filetype."".word(list_L,i)."_steps".steps."_period".period."_iterations".word(list_iter,i)."_J".J."_V".V."_hz".hz."_kick".kick.".txt" every density u 1:(($2)*(-1)**($1/period)) w l title "L =".word(list_L,i)


###### Averages and Fluctuations

reset session

list_L = "4 6 8 10 12"
list_iter = "2560 1280 640 320 160"
#list_iter = "100 100 100 100 100"
list_iter = "320 320 320 320 320"
num_L = words(list_L)

period = "1.00"
J = "0.04"
V = "0.50"
hz = "10.00"
kick = "0.04"
steps = "10000"
param = "period".period."_J".J."_V".V."_hz".hz."_kick".kick
filetype = "data/magnetizations/Sz0_DENSE_SWAP_hz_V_Disorder_Neel_AVG_FLUCT_Imbalance_"
#set title "T = ".period.", J = ".J.", V = ".V.", hz = ".hz.", kick = ".kick
set ylabel "$\\overline{I}$"
set xlabel "$t/T$"
set logscale x
set yrange [-0.1:1.1]
set xrange [1:steps]
set xtics 1, 10, 1e6 format "$10^{%T}$"
set key left bottom box 3
set key at graph 0.04, 0.05

#set terminal png
set size ratio 0.4
set terminal tikz color standalone font ",15"
set output "figures/Scaling_Imbalance_Evolution_Sz0_SWAP_hz_V_Disorder_".param.".tex"
plot for[i=1:num_L] filetype."nspin".word(list_L,i)."_steps".steps."_period".period."_iterations".word(list_iter,i)."_J".J."_V".V."_hz".hz."_kick".kick.".txt" u 1:($2*(-1)**($1/period)) w l lc i title "$L = ".word(list_L,i)."$"


#set logscale y
set yrange [0:0.1]
#set ytics 1e-6, 10, 1e6 format "10^{%T}"
set ytics format "%.2f"
set key left top #box 3
set key at graph 0.04, 0.95

set output "figures/Scaling_Imbalance_Fluctuations_Evolution_Sz0_SWAP_hz_V_Disorder_".param.".tex"
plot for[i=1:num_L] filetype."nspin".word(list_L,i)."_steps".steps."_period".period."_iterations".word(list_iter,i)."_J".J."_V".V."_hz".hz."_kick".kick.".txt" u 1:3 w l lc i title "$L = ".word(list_L,i)."$"
unset logscale y


#Varying J
V = "0.50"
hz = "6.00"
j = 6
#set title "nspin = ".word(list_L,j).", n_{iter} = ".word(list_iter2,j).", period = ".period.", V = ".V.", hz = ".hz.", kick = ".kick
#set output "figures/Avg_Imbalance_Evolution_Sz0_SWAP_hz_Disorder_nspin".word(list_L,j)."_period".period."_V".V."_hz".hz."_kick".kick.".png"

#plot for [i=1:6] filetype."".word(list_L,j)."_steps".steps."_period".period."_iterations".word(list_iter1,j)."_J".word(list_J1,i)."_V".V."_hz".hz."_kick".kick.".txt" every density u 1:(abs($2)) w l title "Swap, J =".word(list_J1,i)

#plot for [i=1:6] filetype."".word(list_L,j)."_steps".steps."_period".period."_iterations".word(list_iter2,j)."_J".word(list_J2,i)."_V".V."_hz".hz."_kick".kick.".txt" every density::::steps u 1:(abs($2)) w l title "Swap, J =".word(list_J2,i)



#Varying V
J = "0.05"
hz = "6.00"
j = 5
#unset logscale x

#set title "nspin = ".word(list_L,j).", n_{iter} = ".word(list_iter3,j).", period = ".period.", J = ".J.", hz = ".hz.", kick = ".kick
#set output "figures/Avg_Imbalance_Evolution_Sz0_SWAP_hz_Disorder_nspin".word(list_L,j)."_period".period."_J".J."_hz".hz."_kick".kick."_varying_V_log.png"

#plot for [i=1:3] filetype."".word(list_L,j)."_steps".steps."_period".period."_iterations".word(list_iter1,j)."_J".J."_V".word(list_V1,i)."_hz".hz."_kick".kick.".txt" every density u 1:(abs($2)) w l title "Swap, V =".word(list_V1,i)

#plot for [i=2:7] filetype."".word(list_L,j)."_steps".steps."_period".period."_iterations".word(list_iter2,j)."_J".J."_V".word(list_V2,i)."_hz".hz."_kick".kick.".txt" every density::::steps  u 1:(abs($2)) w l title "Swap, V =".word(list_V2,i)

#plot for [i=1:20] filetype."".word(list_L,j)."_steps".steps."_period".period."_iterations".word(list_iter3,j)."_J".J."_V".word(list_V3,i)."_hz".hz."_kick".kick.".txt" every density::::steps  u 1:(abs($2)) w l title "Swap, V =".word(list_V3,i)



#Varying hz
J = "0.10"
V = "0.50"
j = 6
#unset logscale x

#set title "nspin = ".word(list_L,j).", n_{iter} = ".word(list_iter2,j).", period = ".period.", J = ".J.", V = ".V.", kick = ".kick
#set output "figures/Avg_Imbalance_Evolution_Sz0_SWAP_hz_Disorder_nspin".word(list_L,j)."_period".period."_J".J."_V".V."_kick".kick."_log.png"

#plot for [i=1:3] filetype."".word(list_L,j)."_steps".steps."_period".period."_iterations".word(list_iter1,j)."_J".J."_V".V."_hz".word(list_hz1,i)."_kick".kick.".txt" every density u 1:(abs($2)) w l title "Swap, hz =".word(list_hz1,i)

#plot for [i=2:7] filetype."".word(list_L,j)."_steps".steps."_period".period."_iterations".word(list_iter2,j)."_J".J."_V".V."_hz".word(list_hz2,i)."_kick".kick.".txt" every density::::steps  u 1:(abs($2)) w l title "Swap, hz =".word(list_hz2,i)


#Varying kick
J = "0.10"
V = "0.50"
hz = "6.00"
j = 5
#unset logscale x

#set title "nspin = ".word(list_L,j).", n_{iter} = ".word(list_iter2,j).", period = ".period.", J = ".J.", V = ".V.", hz = ".hz
#set output "figures/Avg_Imbalance_Evolution_Sz0_SWAP_hz_Disorder_nspin".word(list_L,j)."_period".period."_J".J."_V".V."_hz".hz."_log.png"

#plot for [i=1:5] filetype."".word(list_L,j)."_steps".steps."_period".period."_iterations".word(list_iter2,j)."_J".J."_V".V."_hz".hz."_kick".word(list_kick,i).".txt" every density::::steps  u 1:(abs($2)) w l title "kick =".word(list_kick,i)





### Plot abs(AVG) vs FLUCT 
nspin = '8'
steps = '4000'
n_iter = '30'
kick = '0.00'
Jint = '0.00'
Vint = '1.00'
hx = '0.00'
hz = '3.01'
#set xrange [0.05:2e2]
#set yrange [-0.2:1.1]
#set title 'nspin = '.nspin.', hz = '.hz.', J = '.Jint
#set output "figures/avg_fluct_nspin".nspin."_hz".hz."_J".Jint."_Vint".Vint.".png"
#plot  'filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J'.Jint.'_V'.Vint.'_h'.hx.'_hz'.hz.'_no_kick'.kick.'.txt' u 2:(-$1) w l title 'AVG',\
      'filetype.'_FLUCT_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J'.Jint.'_V'.Vint.'_h'.hx.'_hz'.hz.'_no_kick'.kick.'.txt' u 2:1 w l title 'FLUCT'



### Plot AVG over few steps
nspin = '8'
h_coupling = '0.15'
kick = '0.05'
steps = 100
#set yrange [-1.1:1.1]
#unset logscale x
#set xtics 10 format "%.0s"
#set title 'nspin = '.nspin.', h = '.h_coupling.', eps = '.kick
#set output "figures/avg_steps".steps."_nspin".nspin."_h".h_coupling."_eps".kick.".png"
#plot [1:steps] 'filetype.'_AVG_nspin'.nspin.'_steps100000_iterations40_h'.h_coupling.'_kick'.kick.'.txt' u 2:1 w fillsteps title 'L ='.nspin lt 1




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



### Plot at fixed 'nspin' for varying parameter (J, V, hx, hz)
nspin = '8'
steps = '4000'
n_iter = '30'
kick = '0.00'
Jint = '2.00'
Vint = '2.00'
hx = '0.00'
hz = '0.60'
#set yrange [0:1.1]
#set xrange [0.05:2e2]
set xtics (0.1, 1, 10, 1e2, 1e3, 1e4, 1e5, 1e6) format "10^{%T}"
#set output "figures/avg_imbalance_nspin".nspin."_Vint".Vint."_Dense_MBL_phase.png"
#set output "figures/figure.png"
#set title 'no kick, |1010...> state, nspin = '.nspin.', n_{iter} = '.n_iter.', V = '.Vint
set ylabel 'I(t)'
set xlabel 't(in units of hbar/2J)'
#plot for [hz in "1.00 2.00 3.00 4.00 6.00"] 'filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J'.Jint.'_V'.Vint.'_h'.hx.'_hz'.hz.'_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz ='.hz

#plot  'filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J1.00_V'.Vint.'_h'.hx.'_hz3.00_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = 3, J = 1',\
      'filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J0.50_V'.Vint.'_h'.hx.'_hz3.00_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = 3, J = 0.5',\
      'filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J1.00_V'.Vint.'_h'.hx.'_hz4.00_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = 4, J = 1',\
      'filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J2.00_V'.Vint.'_h'.hx.'_hz6.00_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = 6, J = 2'

#set xrange [0.05:2e1]
#set yrange [-0.2:1.1]
#set output "figures/avg_imbalance_nspin".nspin."_Vint".Vint."_Dense_Ergodic_phase.png"
#plot  'filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J1.50_V'.Vint.'_h'.hx.'_hz3.00_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = 3, J = 1.5',\
      'filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J2.00_V'.Vint.'_h'.hx.'_hz0.00_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = 0, J = 2'

Vint = '1.00'
Jint = '0.00'
#set title 'no kick, |1010...> state, nspin = '.nspin.', n_{iter} = '.n_iter.', V = '.Vint.', Jint = '.Jint
#set output "figures/avg_imbalance_nspin".nspin."_Vint".Vint."_Dense_AL_phase.png"
#plot for [ hz in "2.50 3.00 4.00" ] 'filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J'.Jint.'_V'.Vint.'_h'.hx.'_hz'.hz.'_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = '.hz


#filename = "data/magnetizations/Sz0_vs_Full_MBL_hz_Disorder_SPARSE_AVG_FLUCT_Imbalance_nspin10_steps100000_iterations1_J0.50_V1.00_hz3.00_kdim30.txt"
#filename = "data/magnetizations/Sz0_vs_Full_MBL_hz_Disorder_SPARSE_AVG_FLUCT_Imbalance_nspin10_steps4000_iterations30_J0.50_V1.00_hz3.00_kdim30.txt"
filename = "data/magnetizations/Sz0_SPARSE_vs_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin12_steps12000_iterations1_J0.50_V1.00_hz3.00_kdim30.txt"


reset session
