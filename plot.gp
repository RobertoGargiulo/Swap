set terminal png # tikz color # standalone
set output "figures/figure.png"
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
title2 = 'L = 2'
title3 = 'L = 3'
title4 = 'L = 4'
title5 = 'L = 5'
type = 'AVG'


### Plot abs(AVG) or FLUCT at varying 'nspin'
plot  'data/magnetizations/Swap_Sz_'.type.'_nspin2_steps100000_iterations40_h0.30_kick0.10.txt' u 2:(abs($1)) w l title title2 lt 1 ,\
      'data/magnetizations/Swap_Sz_'.type.'_nspin3_steps100000_iterations40_h0.30_kick0.10.txt' u 2:(abs($1)) w l title title3 lt 2 ,\
      'data/magnetizations/Swap_Sz_'.type.'_nspin4_steps100000_iterations40_h0.30_kick0.10.txt' u 2:(abs($1)) w l title title4 lt 3 ,\
      'data/magnetizations/Swap_Sz_'.type.'_nspin5_steps1000000_iterations40_h0.30_kick0.10.txt' u 2:(abs($1)) w l title title5 lt 4


nspin = '5'
### Plot abs(AVG) vs FLUCT 

# plot  'data/magnetizations/Swap_Sz_AVG_nspin'.nspin.'_steps1000000_iterations40_h0.30_kick0.10.txt' u 2:(abs($1)) w l title 'AVG' lt 1 ,\
#      'data/magnetizations/Swap_Sz_FLUCT_nspin'.nspin.'_steps1000000_iterations40_h0.30_kick0.10.txt' u 2:(abs($1)) w l title 'FLUCT' lt 2

### Plot AVG
#plot 'data/magnetizations/Swap_Sz_AVG_nspin2_steps100000_iterations40.txt' u 2:1 w fillsteps title title2 lt 1
#plot 'data/magnetizations/Swap_Sz_AVG_nspin5_steps1000000_iterations40.txt' u 2:1 w fillsteps title title5 lt 2


### Plot at fixed 'nspin' for varying 'h' or 'eps'
# plot for [i in "0.00 0.10 0.40 0.70 1.00"] 'data/magnetizations/Swap_Sz_AVG_nspin4_steps100000_iterations40_h0.30_kick'.i.'.txt' u 2:(abs($1)) w l title 'eps ='.i


set output
unset logscale x
