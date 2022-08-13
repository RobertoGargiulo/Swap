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

#list(start,end,increment)=system(sprintf("seq %g %g %g", start, increment, end))
start=0.; end=6.; inc=0.6
list = ""; a=start-inc; while (a<end) {list=list.sprintf(" %.2f",a=a+inc)}
list = "6 8 10 12 14"

set yrange [-0.1:1.1]
set ytics format "%.1f"
set xtics format "%.1f"
set xlabel 'h_{z}/J'
set title 'no kick, |1010...> state, n_{iter} = '.iter.', J = '.J.' V = '.V

#Imbalance

#set ylabel 'Imbalance'


##########

file_close2 = "data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin"
#set output "figures/MBL_PT_Scaling_Imbalance_J".J."_V".V."_Dense_close_to_hz2.png"
#plot for [nspin in list] file_close2."".nspin."_time_step".time_step."_steps".steps."_close_to_hz2.txt" u ($3/$1):(-$4):5 w lp lw 1.6 title 'L ='.nspin

#plot for [nspin in list] file_close2."".nspin."_time_step".time_step."_steps".steps."_close_to_hz2.txt" u ($3/$1):(-$4):5 w errorlines title 'L ='.nspin

##########

file_upto6 = "data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin"

set xrange [1:12]
set logscale x
set xtics 0.5, 2, 16 format "%.2f"

list1 = "6 8 10 12 14 16"
list2 = "640 320 160 80 40 20"
set title "no kick, |1010...> state, J = ".J." V = ".V
#set output "figures/MBL_PT_Scaling_Imbalance_J".J."_V".V."_Dense_up_to_hz6.png"
#plot for [i=1:5] file_upto6."".word(list1,i)."_time_step".time_step."_steps".steps."_iterations".word(list2,i)."_J".J."_V".V."_up_to_hz6.txt" u ($3/$1):(-$4):5 w lp lw 1.6 title "L = ".word(list1,i)."  n_{iter} = ".word(list2,i)

#set output "figures/MBL_PT_Scaling_Imbalance_J".J."_V".V."_Dense_up_to_hz6_w_errors.png"
#plot for [i=1:4] file_upto6."".word(list1,i)."_time_step".time_step."_steps".steps."_iterations".word(list2,i)."_J".J."_V".V."_up_to_hz6.txt" u ($3/$1):(-$4):5 w errorlines title "L = ".word(list1,i)."  n_{iter} = ".word(list2,i)


###############

set xrange [0.2:40]
file_imb_upto10="data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin"

#set output "figures/MBL_PT_Scaling_Imbalance_J".J."_V".V."_Dense_up_to_hz10_w_errors.png"
#plot for [i=1:4] file_imb_upto10."".word(list1,i)."_time_step".time_step."_steps".steps."_iterations".word(list2,i)."_J".J."_V".V."_up_to_hz10.txt" u ($3/$1):(-$4):5 w errorlines title "L = ".word(list1,i)."  n_{iter} = ".word(list2,i)








#Average Gap Ratio

set yrange [0.35:0.6]
set ytics format "%.2f"
set ylabel 'Average Gap Ratio'

set xrange [0.2:12]
set logscale x
set xtics 0.5, 2, 8 format "%.2f"

file_close2 = "data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin"
#set output "figures/MBL_PT_Scaling_Gap_Ratio_J".J."_V".V."_Dense_close_to_hz2.png"
#plot for [nspin in list] file_close2."".nspin."_time_step".time_step."_steps".steps."_close_to_hz2.txt" u ($3/$1):6 w lp lw 1.6 title 'L ='.nspin

#set output "figures/MBL_PT_Scaling_Gap_Ratio_J".J."_V".V."_Dense_close_to_hz2_w_errors.png"
#plot for [nspin in list] file_close2."".nspin."_time_step".time_step."_steps".steps."_close_to_hz2.txt" u ($3/$1):6:7 w errorlines title 'L ='.nspin



#file_upto10 = "data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_Gap_Ratio_nspin"
#set output "figures/MBL_PT_Scaling_Gap_Ratio_J".J."_V".V."_Dense_up_to_hz10.png"
#plot for [nspin in list] file_upto10."".nspin."_iterations80_J".J."_V".V."_up_to_hz10.txt" u ($3/$1):($4) w l title 'L ='.nspin


###############

file_upto6 = "data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin"

set xrange [1:12]
set logscale x
set xtics 1, 2, 8 format "%.2f"

list1 = "6 8 10 12 14 16"
list2 = "640 320 160 80 40 20"
set title "no kick, J = ".J." V = ".V

#set output "figures/MBL_PT_Scaling_Gap_Ratio_J".J."_V".V."_Dense_up_to_hz6.png"
#plot for [i=1:5] file_upto6."".word(list1,i)."_time_step".time_step."_steps".steps."_iterations".word(list2,i)."_J".J."_V".V."_up_to_hz6.txt" u ($3/$1):6 w lp lw 1.6 title "L = ".word(list1,i)."  n_{iter} = ".word(list2,i)

##############


#Fixed Disorder, Varying 'nspin'

start=0.00
inc=0.10
end=1.00
unset xrange
unset logscale x
set xtics 6, 2, 14
set xlabel "L"
set ylabel "Imbalance"
#list = ""; a=start-inc; while (a<end) {list=list.sprintf(" %.2f",a=a+inc)}
#system(sprintf(list))

J="0.50"
V="0.25"
set title "no kick, |1010...> state, J = ".J." V = ".V
file_upto_nspin14="data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_time_step0.50_steps4000_J0.50_V0.25_hz"
set yrange [-0.1:1.1]
#set output "figures/MBL_PT_Scaling_Imbalance_J".J."_V".V."_Dense_up_to_nspin14.png"
#plot for [hz in "0.05 1.00 3.50 6.00 10.00 30.00"] file_upto_nspin14."".hz."_up_to_nspin12.txt" u 1:(-$4):5 w errorlines lw 1.6 title "hz = ".hz

#Fixed nspin,J,hz,kick,iter,steps, varying V

file="Swap_varying_V_data2.txt"
file="Swap_varying_V.txt"

set ylabel "Imbalance"
set xlabel "V"
set xtics 0, pi/2, 4*pi

#mod(x) = system( sprintf("echo $((%d % 4))",x)  )
#system(sprintf("echo %s", mod(4)))
#do for [i=1:10] {set xtics add (sprintf("%d/4*pi",i) i)}
set yrange [0:1.1]

nspin = "10"
iter = "160"
period = "1.00"
J = "0.05"
hz = "6.00"
kick = "0.00"

set title "nspin = ".nspin.", n_{iter} = ".iter.", period = ".period.", J = ".J.", hz = ".hz.", kick = ".kick
set output "figures/Avg_Imbalance_Sz0_SWAP_hz_Disorder_nspin".nspin."_period".period."_J".J."_hz".hz."_kick".kick."_up_to_V12.png"
plot file u 3:(-$7) w lp lw 1.4 title "Swap", "" u 3:(-$9) title "MBL"


##### Decay Time


#Fixed nspin,J,hz,kick,iter,steps, varying V

file_time="Swap_decay_times.txt"

set ylabel "Decay Time"
set xlabel "L"
set logscale y
unset logscale x


list_L = "2 4 6 8 10"
list_iter = "5120 2560 1280 640 320"
list_J = "0.02 0.04 0.06 0.08 0.10"
list_V = "0.25 0.50"
list_h = "3.00 6.00 12.00"

period = "1.00"
kick = "0.00"

set title "n_{iter} = 2^{1-L/2}5120, steps = 10^{5}, period = ".period.", J = ".J.", hz = ".hz.", kick = ".kick
set output "figures/Decay_Time_Sz0_SWAP_hz_Disorder_nspin".nspin."_period".period."_J".J."_V".V."_kick".kick."_up_to_L10.png"
plot file u 1:11 every 3:20:1 w lp lw 1.4 title "hz = ".hz #, "" u 3:(-$9) title "MBL"





set output
