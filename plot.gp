set terminal png # tikz color # standalone
#set output "figures/figure.png"
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
# set label 1 at graph 0.5, 0.75 "$k$" center


# filetype = 'Swap_Sz'
filetype = 'Clean_MBL_Imbalance'


### Plot abs(AVG) or FLUCT at varying 'nspin'
title1 = 'L = 2'
title2 = 'L = 4'
title3 = 'L = 6'
title4 = 'L = 8'
type = 'AVG'
h_coupling = '0.15'
kick = '0.05'
#set output "figures/".type."_h".h_coupling."_eps".kick.".png"
#set title "mag ".type.", h = ".h_coupling.", eps = ".kick
#plot  'data/magnetizations/'.filetype.'_'.type.'_nspin2_steps100000_iterations40_h'.h_coupling.'_kick'.kick.'.txt' u 2:(abs($1)) w l title title1 lt 1 ,\
#      'data/magnetizations/'.filetype.'_'.type.'_nspin4_steps100000_iterations40_h'.h_coupling.'_kick'.kick.'.txt' u 2:(abs($1)) w l title title2 lt 2 ,\
#      'data/magnetizations/'.filetype.'_Sz_'.type.'_nspin6_steps100000_iterations40_h'.h_coupling.'_kick'.kick.'.txt' u 2:(abs($1)) w l title title3 lt 3 ,\
#      'data/magnetizations/'.filetype.'_'.type.'_nspin8_steps100000_iterations40_h'.h_coupling.'_kick'.kick.'.txt' u 2:(abs($1)) w l title title4 lt 4 ,\



### Plot abs(AVG) vs FLUCT 
nspin = '8'
steps = '4000'
n_iter = '30'
kick = '0.00'
Jint = '0.00'
Vint = '1.00'
hx = '0.00'
hz = '3.01'
set xrange [0.05:2e2]
set yrange [-0.2:1.1]
set title 'nspin = '.nspin.', hz = '.hz.', J = '.Jint
set output "figures/avg_fluct_nspin".nspin."_hz".hz."_J".Jint."_Vint".Vint.".png"
plot  'data/magnetizations/'.filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J'.Jint.'_V'.Vint.'_h'.hx.'_hz'.hz.'_no_kick'.kick.'.txt' u 2:(-$1) w l title 'AVG',\
      'data/magnetizations/'.filetype.'_FLUCT_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J'.Jint.'_V'.Vint.'_h'.hx.'_hz'.hz.'_no_kick'.kick.'.txt' u 2:1 w l title 'FLUCT'



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
#plot [1:steps] 'data/magnetizations/'.filetype.'_AVG_nspin'.nspin.'_steps100000_iterations40_h'.h_coupling.'_kick'.kick.'.txt' u 2:1 w fillsteps title 'L ='.nspin lt 1




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
#plot [1:steps] 'data/magnetizations/'.filetype.'_nspin'.nspin.'_steps100000_iterations40_h0.30_kick0.10.txt' every ::start::start+steps u ($2-start):1 w fillsteps title 'iteration = 1',\
#              'data/magnetizations/'.filetype.'_nspin'.nspin.'_steps100000_iterations40_h0.30_kick0.10.txt' every ::1e5+1::1e5+steps u 2:1 w fillsteps title 'iteration = 2'



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
set title 'no kick, |1010...> state, nspin = '.nspin.', n_{iter} = '.n_iter.', V = '.Vint
set ylabel 'I(t)'
set xlabel 't(in units of hbar/V)'
#plot for [hz in "1.00 2.00 3.00 4.00 6.00"] 'data/magnetizations/'.filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J'.Jint.'_V'.Vint.'_h'.hx.'_hz'.hz.'_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz ='.hz

#plot  'data/magnetizations/'.filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J1.00_V'.Vint.'_h'.hx.'_hz3.00_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = 3, J = 1',\
      'data/magnetizations/'.filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J0.50_V'.Vint.'_h'.hx.'_hz3.00_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = 3, J = 0.5',\
      'data/magnetizations/'.filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J1.00_V'.Vint.'_h'.hx.'_hz4.00_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = 4, J = 1',\
      'data/magnetizations/'.filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J2.00_V'.Vint.'_h'.hx.'_hz6.00_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = 6, J = 2'

#set xrange [0.05:2e1]
#set yrange [-0.2:1.1]
#set output "figures/avg_imbalance_nspin".nspin."_Vint".Vint."_Dense_Ergodic_phase.png"
#plot  'data/magnetizations/'.filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J1.50_V'.Vint.'_h'.hx.'_hz3.00_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = 3, J = 1.5',\
      'data/magnetizations/'.filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J2.00_V'.Vint.'_h'.hx.'_hz0.00_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = 0, J = 2'

Vint = '1.00'
Jint = '0.00'
#set title 'no kick, |1010...> state, nspin = '.nspin.', n_{iter} = '.n_iter.', V = '.Vint.', Jint = '.Jint
#set output "figures/avg_imbalance_nspin".nspin."_Vint".Vint."_Dense_AL_phase.png"
#plot for [ hz in "2.50 3.00 4.00" ] 'data/magnetizations/'.filetype.'_AVG_nspin'.nspin.'_steps'.steps.'_iterations'.n_iter.'_J'.Jint.'_V'.Vint.'_h'.hx.'_hz'.hz.'_no_kick'.kick.'.txt' u 2:(-$1) w l title 'hz = '.hz


set output
unset logscale x
