set terminal png # tikz color # standalone
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

############## Varying J

file_swap_J="sort_Swap_entanglement.txt"

unset logscale x
list_hz = "0.01 2.00 4.00 8.00 12.00 16.00"
list_J = "0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50"
V = "0.50"
hz = "6.00"
kick = "0.00"
period = "1.00"

list_L = "4 6 8 10"
num_L = words(list_L)
num_J = words(list_J)

set ylabel "Gap Ratio"
set title "V = ".V.", hz = ".hz." kick = ".kick.", period = ".period


set xlabel "J"
set xtics 0, 0.10, 0.50 format "%.2f"
set xrange [0:0.50]
set yrange [0.35:0.6]
#set output "figures/Swap_PT_Scaling_Gap_Ratio_V".V."_hz".hz."_up_to_J0.5.png"
#plot for [i=0:num_L-1] file_swap_J u 2:7:8 every :::i::i w errorlines title "L = ".word(list_L,i+1)

set xlabel "L"
set xtics 4, 2, 10 format "%.0f"
set xrange [4:10]
set yrange [0.35:0.6]
#set output "figures/Swap_PT_Scaling_Gap_Ratio_V".V."_hz".hz."_up_to_J0.5.png"
#plot for [i=0:(num_J-1)/2] file_swap_J u 1:7:8 every :::2*i::2*i w errorlines title "J = ".word(list_J,i+1)

set xlabel "L"
set xtics 4, 2, 10 format "%.0f"
set xrange [4:10]
unset yrange
#set output "figures/Swap_PT_Log_Gap_V".V."_hz".hz."_up_to_J0.5.png"
#plot for [i=0:(num_J-1)/2] file_swap_J u 1:11:12 every :::2*i::2*i w errorlines title "J = ".word(list_J,i+1)

set xlabel "J"
set xtics 0, 0.1, 0.5 format "%.2f"
set xrange [0:0.5]

#file_swap_J='< ./spaces.sh Swap_entanglement.txt'
#set output "figures/Swap_PT_Log_Near_Pair_Gap_V".V."_hz".hz."_up_to_J0.5.png"
#plot for [i=0:num_L-1] file_swap_J u 2:13 every :::i::i w errorlines pt 4 lt i+1 title "log Delta, L = ".word(list_L,i+1),\
#     for [i=0:num_L-1] file_swap_J u 2:14 every :::i::i w errorlines pt 5 lt i+1 notitle #"log Delta_0, L = ".word(list_L,i+1)


list_L = "4 6 8 10 12"
num_L = words(list_L)
set xtics (0.01, 0.05, 0.1, 0.5, 1) format "%.2f"
set xrange [0.01:1]
file_swap_J="sort_Swap_spectrum_J.txt"
set logscale x
set key right bottom box 3
#set output "figures/Swap_PT_Log_Near_Pair_Gap_V".V."_hz".hz."_up_to_J1.png"
#plot for [i=0:num_L-1] file_swap_J u 2:13 every :::i::i w errorlines pt 4 lt i+1 title "log Delta, L = ".word(list_L,i+1),\
#     for [i=0:num_L-1] file_swap_J u 2:14 every :::i::i w errorlines pt 5 lt i+1 notitle #"log Delta_0, L = ".word(list_L,i+1)


set xlabel "h_{z}"
set xtics format "%.2f"
set xrange [0.01:16.00]
file_swap_hz="sort_Swap_spectrum_hz.txt"
set logscale x
set key left bottom box 3
#set output "figures/Swap_PT_Log_Near_Pair_Gap_J".J."V".V."_up_to_hz16.png"
#plot for [i=0:num_L-1] file_swap_hz u 4:13 every :::i::i w errorlines pt 4 lt i+1 title "log Delta, L = ".word(list_L,i+1),\
#     for [i=0:num_L-1] file_swap_hz u 4:14 every :::i::i w errorlines pt 5 lt i+1 notitle #"log Delta_0, L = ".word(list_L,i+1)




########


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
#set output "figures/Avg_Imbalance_Sz0_SWAP_hz_Disorder_nspin".nspin."_period".period."_J".J."_hz".hz."_kick".kick."_up_to_V12.png"
#plot file u 3:(-$7) w lp lw 1.4 title "Swap", "" u 3:(-$9) title "MBL"


##### Decay Time


#Fixed nspin,J,hz,kick,iter,steps, varying V

file_time="sort_Swap_decay_times2.txt"
file_time_J="sort_Swap_decay_times_varying_J.txt"
file_time_V="sort_Swap_decay_times_varying_V.txt"
file_time_hz="sort_Swap_decay_times_varying_hz.txt"
file_time_eps="sort_Swap_decay_times_varying_eps.txt"
#file_time2="Swap_decay_times2.txt"


set ylabel "Decay Time"
set xlabel "L"
set logscale y
unset logscale x


list_L = "2 4 6 8 10 12"
list_iter = "5120 2560 1280 640 320 160"
list_J = "0.02 0.04 0.06 0.08 0.10"
list_V = "0.25 0.50"
list_hz = "3.00 6.00 12.00"

list_J = "0.00 0.04 0.08 0.12 0.16 0.20 0.24"
hz = "6.00"
V = "0.25"



period = "1.00"
kick = "0.00"
steps = 1e6
set xrange [1.50:10.50]
set yrange [1:steps]

set xtics list_L
set ytics 0.001, 10, 10e10 format "10^{%T}"

# J = 0.02, ..., 0.10, V = 0.25, hz = 3.00 -> j = 0, ..., 4
# J = 0.02, ..., 0.10, V = 0.50, hz = 3.00 -> j = 5, ..., 9
# J = 0.02, ..., 0.10, V = 0.25, hz = 6.00 -> j = 10, ..., 14
# J = 0.02, ..., 0.10, V = 0.50, hz = 6.00 -> j = 15, ..., 19
# J = 0.02, ..., 0.10, V = 0.25, hz = 12.00 -> j = 20, ..., 24
# J = 0.02, ..., 0.10, V = 0.50, hz = 12.00 -> j = 25, ..., 29
 
# INDEX = i+5*j+10*k -> 0<i<4 (changing J), 0<j<1 (changing V), 0<k<2 (changing hz)
# In generale: i + num_J * j + num_J * num_V * k, 0<i<num_J-1; 0<j<num_V-1; 0<k<num_hz-1


j = 0; k = 0
#set title "n_{iter} = 2^{1-L/2}5120, steps = 10^{6}, period = ".period.", V = ".word(list_V,j+1).", hz = ".word(list_hz,k+1).", kick = ".kick

list_J = "0.05 0.10 0.15 0.20 0.25 0.30"
num_J = words(list_J)
num_V = 1
num_hz = 1
V = "0.25"
hz = "6.00"
set yrange [1:1e9]
set xrange [2:12]
set key left top box 3
set title "n_{iter} = 2^{1-L/2}5120, period = ".period.", V = ".V.", hz = ".hz.", kick = ".kick
#set output "figures/Scaling_Decay_Time_varying_J.png"
#plot for [i=0:num_J-1] file_time_J every :::i+num_J*j+num_J*num_V*k::i+num_J*j+num_J*num_V*k u 1:7:8 w errorlines title "J = ".word(list_J,i+1)

i = 0; k = 0
list_V = "0.00 0.40 0.80 1.20 1.60"
J = "0.10"
hz = "6.00"
kick = "0.00"
num_V = words(list_V)
num_J = 1
num_hz = 1
#set title "n_{iter} = 2^{1-L/2}5120, steps = 10^{5}, period = ".period.", J = ".word(list_J,i+1).", hz = ".word(list_hz,k+1).", kick = ".kick
set title "n_{iter} = 2^{1-L/2}5120, steps = 10^{6}, period = ".period.", J = ".J.", hz = ".hz.", kick = ".kick
#set output "figures/Scaling_Decay_Time_varying_V.png"
#plot for [j=0:num_V-1] file_time_V every :::i+num_J*j+num_J*num_V*k::i+num_J*j+num_J*num_V*k u 1:9:10 w errorlines title "V = ".word(list_V,j+1)

i = 0; j = 0
list_hz = "0.01 2.00 4.00 8.00 12.00 16.00"
J = "0.10"
V = "0.50"
kick = "0.00"
num_hz = words(list_hz)
num_J = 1
num_V = 1
#set title "n_{iter} = 2^{1-L/2}5120, steps = 10^{5}, period = ".period.", J = ".word(list_J,i+1).", V = ".word(list_V,j+1).", kick = ".kick
set title "n_{iter} = 2^{1-L/2}5120, steps = 10^{6}, period = ".period.", J = ".J.", V = ".V.", kick = ".kick
#set output "figures/Scaling_Decay_Time_varying_hz.png"
#plot for [k=0:num_hz-1] file_time_hz every :::i+num_J*j+num_J*num_V*k::i+num_J*j+num_J*num_V*k u 1:9:10 w errorlines title "hz = ".word(list_hz,k+1)

set yrange [1:1e9]
set xrange [2:12]
set key left top box 3
file_hz="sort_Swap_decay_times_varying_hz.txt"
#set output "figures/Scaling_Decay_Time_varying_hz.png"
#plot for [k=0:num_hz-1] file_hz every :::k::k u 1:7:8 w errorlines title "hz = ".word(list_hz,k+1)

i = 1; j = 0; k = 0
do for [i=0:num_J-1] {
#set table $DataSelected
#  plot file_time_V every :::i+num_J*j+num_J*num_V*k::i+num_J*j+num_J*num_V*k u 1:2:3:4:9:10:11 w table
#unset table
#print "Table"
#print $DataSelected
}



period = "1.00"
V = "0.50"
hz = "6.00"
kick = "0.00"

list_L = "2 4 6 8 10 12 14 16"
list_iter = "2560 1280 640 320 160 80 40 20 10"
list_J = "0.05 0.10 0.15 0.20 0.25 0.30 0.40 0.50 1.00"

init_state = "Neel"
file_time_J = "sort_Swap_".init_state."_decay_times_varying_J.txt"
num_J = words(list_J)
num_V = 1
num_hz = 1
V = "0.25"
hz = "6.00"
set yrange [1:1e9]
set xrange [2-0.3:10+0.3]
set key left top box 3
set title "n_{iter} = 2^{1-L/2} 2560, period = ".period.", J = ".J.", V = ".V.", hz = ".hz.", kick = ".kick."\n{/*0.8 ".init_state." initial state}"
set output "figures/Scaling_Decay_Time_".init_state."_varying_J.png"
plot for [i=0:num_J-1] file_time_J every :::i+num_J*j+num_J*num_V*k::i+num_J*j+num_J*num_V*k u 1:7:8 w errorlines title "J = ".word(list_J,i+1)

set output
