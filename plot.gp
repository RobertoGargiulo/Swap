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
set xrange [1:1e6]
# set key notitle invert under reverse Left left spacing 2 samplen 0.7
# set arrow 1 filled from graph 0.4, 0.7 to graph 0.6, 0.7
# set label 1 at graph 0.5, 0.75 "$k$" center



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
#plot  'data/magnetizations/Swap_Sz_'.type.'_nspin2_steps100000_iterations40_h'.h_coupling.'_kick'.kick.'.txt' u 2:(abs($1)) w l title title1 lt 1 ,\
#      'data/magnetizations/Swap_Sz_'.type.'_nspin4_steps100000_iterations40_h'.h_coupling.'_kick'.kick.'.txt' u 2:(abs($1)) w l title title2 lt 2 ,\
#      'data/magnetizations/Swap_Sz_'.type.'_nspin6_steps100000_iterations40_h'.h_coupling.'_kick'.kick.'.txt' u 2:(abs($1)) w l title title3 lt 3 ,\
#      'data/magnetizations/Swap_Sz_'.type.'_nspin8_steps100000_iterations40_h'.h_coupling.'_kick'.kick.'.txt' u 2:(abs($1)) w l title title4 lt 4 ,\



### Plot abs(AVG) vs FLUCT 
nspin = '2'
h_coupling = '0.30'
kick = '0.10'
#set title 'nspin = '.nspin.', h = '.h_coupling.', eps = '.kick
#set output "figures/avg_fluct_nspin".nspin."_h".h_coupling."_eps".kick.".png"
#plot  'data/magnetizations/Swap_Sz_AVG_nspin'.nspin.'_steps100000_iterations40_h'.h_coupling.'_kick'.kick.'.txt' u 2:(abs($1)) w l title 'AVG' lt 1 ,\
     'data/magnetizations/Swap_Sz_FLUCT_nspin'.nspin.'_steps100000_iterations40_h0.30_kick0.10.txt' u 2:(abs($1)) w l title 'FLUCT' lt 2



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
#plot [1:steps] 'data/magnetizations/Swap_Sz_AVG_nspin'.nspin.'_steps100000_iterations40_h'.h_coupling.'_kick'.kick.'.txt' u 2:1 w fillsteps title 'L ='.nspin lt 1




### Plot single iteration
nspin = '2'
h_coupling = '0.15'
kick = '0.05'
steps = 100
type = 'early'
set yrange [-1.1:1.1]
unset logscale x
set xtics 10 format "%0.s"
set title 'nspin = '.nspin.', h = '.h_coupling.', eps = '.kick
set output "figures/".type."_steps".steps."_nspin".nspin."_h".h_coupling."_eps".kick.".png"
plot [1:steps] 'data/magnetizations/Swap_Sz_nspin'.nspin.'_steps100000_iterations40_h0.30_kick0.10.txt' every ::1::steps u 2:1 w fillsteps title 'iteration = 1',\
#              'data/magnetizations/Swap_Sz_nspin'.nspin.'_steps100000_iterations40_h0.30_kick0.10.txt' every ::1e5+1::1e5+steps u 2:1 w fillsteps title 'iteration = 2'



### Plot at fixed 'nspin' for varying 'h' or 'eps'
nspin = '4'
h_coupling = '0.0'
kick = '0.00'
#set title 'nspin = '.nspin.', h = '.h_coupling
#set title 'nspin = '.nspin.', eps = '.kick
#plot for [i in "0.00 0.05 0.15 0.30 0.60"] 'data/magnetizations/Swap_Sz_AVG_nspin'.nspin.'_steps100000_iterations40_h'.i.'_kick'.kick.'.txt' u 2:(abs($1)) w l title 'h ='.i


set output
unset logscale x
