set terminal png # tikz color # standalone
#set output "figures/figure.png"
set size ratio 0.5
set margins -1,2,0,0

#Scaling Plot

steps="4000"
iter="100"
J = "0.50"
V = "0.25"
time_step="0.50"
filetype = "data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin"

#list(start,end,increment)=system(sprintf("seq %g %g %g", start, increment, end))
start=0.; end=6.; inc=0.6
list = ""; a=start-inc; while (a<end) {list=list.sprintf(" %.2f",a=a+inc)}
list = "6 8 10 12 14"

set yrange [-0.1:1.1]
set ytics format "%.1f"
set xtics format "%.1f"
set xlabel 'h_{z}/J'
set title 'no kick, |1010...> state, n_{iter} = '.iter.', J = '.J.' V = '.V
#set output "figures/MBL_PT_Scaling_Imbalance_J".J."_V".V."_Dense_close_to_hz2.png"
#set ylabel 'Imbalance'
#plot for [nspin in list] filetype."".nspin."_time_step".time_step."_steps".steps."_close_to_hz2.txt" u ($3/$1):(-$4):5 w errorlines title 'L ='.nspin


#filetype = "data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_Gap_Ratio_nspin"
set yrange [0.35:0.6]
#set xrange [1:20]
set ytics format "%.2f"
set ylabel 'Average Gap Ratio'
#set output "figures/MBL_PT_Scaling_Gap_Ratio_J".J."_V".V."_Dense_close_to_hz2.png"
#plot for [nspin in list] filetype."".nspin."_time_step".time_step."_steps".steps."_close_to_hz2.txt" u ($3/$1):6 w l title 'L ='.nspin
set output "figures/MBL_PT_Scaling_Gap_Ratio_J".J."_V".V."_Dense_close_to_hz2_w_errors.png"
plot for [nspin in list] filetype."".nspin."_time_step".time_step."_steps".steps."_close_to_hz2.txt" u ($3/$1):6:7 w errorlines title 'L ='.nspin
#set output "figures/MBL_PT_Scaling_Gap_Ratio_J".J."_V".V."_Dense_up_to_hz10.png"
#plot for [nspin in list] filetype."".nspin."_iterations80_J".J."_V".V."_up_to_hz10.txt" u ($3/$1):($4) w l title 'L ='.nspin



set output
