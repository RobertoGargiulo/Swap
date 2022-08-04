set terminal png # tikz color # standalone
#set output "figures/figure.png"
set size ratio 0.5
set margins -1,2,0,0

#Scaling Plot

steps="4000"
iter="80"
J = "0.50"
V = "0.25"
time_step="0.50"
filetype = "data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin"
#list(start,end,increment)=system(sprintf("seq %g %g %g", start, increment, end))
start=0.; end=6.; inc=0.6
list = ""; a=start-inc; while (a<end) {list=list.sprintf(" %.2f",a=a+inc)}
list = "6 8 10 12" # 14

set yrange [-0.1:1.1]
set ytics format "%.1f"
set xtics format "%.1f"
set output "figures/MBL_PT_Scaling_Imbalance_J".J."_V".V."_Dense_close_to_hz1.50.png"
set title 'no kick, |1010...> state, n_{iter} = '.iter.', J = '.J.' V = '.V
set ylabel 'Imbalance'
set xlabel 'h_{z}/J'
#plot for [nspin in list] filetype."".nspin."_time_step".time_step."_steps".steps."_close_to_hz1.50.txt" u ($3/$1):(-$4) w l title 'L ='.nspin

plot for [nspin in list] filetype."".nspin."_time_step".time_step."_steps".steps."_close_to_hz1.50.txt" u ($3/$1):($6) w l title 'L ='.nspin



set output
