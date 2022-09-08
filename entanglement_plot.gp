set terminal png # tikz color # standalone
set size ratio 0.5
set margins -1,2,0,0

max(x, y) = (x > y ? x : y)
min(x, y) = (x < y ? x : y)


###### Parameters
nspin = "12"
n_lim = min(nspin/2, 2)
nspin_A = 2
period="1.00"
J = "0.00"
V = "0.50"
hz = "16.00"
kick = "0.00"
######
#unset key

file="data/eigen/Sz0_DENSE_SWAP_hz_Disorder_Entanglement_nspin".nspin."_period".period."_iterations1_J".J."_V".V."_hz".hz."_kick".kick.".txt"
file1="data/eigen/Sz0_DENSE_SWAP_hz_Disorder_Entanglement_nspin".nspin."_period".period."_iterations1_J0.00_V".V."_hz".hz."_kick".kick.".txt"
file2="data/eigen/Sz0_DENSE_SWAP_hz_Disorder_Entanglement_nspin".nspin."_period".period."_iterations1_J0.05_V".V."_hz".hz."_kick".kick.".txt"
figure_param = "_nspin".nspin."_nspinA".nspin_A."_period".period."_iterations1_J".J."_V".V."_hz".hz."_kick".kick.".png"
figure_param1 = "_nspin".nspin."_nspinA".nspin_A."_period".period."_iterations1_J0.00_V".V."_hz".hz."_kick".kick.".png"
figure_param2 = "_nspin".nspin."_nspinA".nspin_A."_period".period."_iterations1_J0.05_V".V."_hz".hz."_kick".kick.".png"

set title "nspin = ".nspin.", J = ".J.", V = ".V.", hz = ".hz.", kick = ".kick.", period = ".period
set yrange [-0.1:nspin_A*log(2)+0.1]
set xlabel "alpha" 
set ylabel "Mutual Information, L_A = ".nspin_A
#set output "figures/Swap_Entanglement_MI".figure_param
#plot file u nspin_A




################ Only LI vs MI (at different J)
set yrange [-0.1:2.1]
set xrange [-0.1:nspin_A*log(2)+0.1]
set xlabel "Mutual Information L_A =".nspin_A 
set ylabel "Local Imbalance" 
#set output "figures/figure.png"
#plot file1 u nspin_A:n_lim+1 title "J = 0.00", file2 u nspin_A:n_lim+1 title "J = 0.05"
################


##### LI, IPR, MBE, I^2 vs MI

plot_title = "nspin = ".nspin.", J = ".J.", V = ".V.", hz = ".hz.", kick = ".kick.", period = ".period

#set output "figures/Swap_Entanglement_MI_corr".figure_param
unset title
set xtics 0, nspin_A*log(2)/3, nspin_A*log(2) format "%.1f"

#set multiplot layout 2,2 columnsfirst scale 1,1 title plot_title

set xrange [-0.1:nspin_A*log(2)+0.1]
set yrange [-0.1:2.1]
set xlabel "Mutual Information, L_A = ".nspin_A
set ylabel "Local Imbalance"
#"I^" #"MBE" #"IPR" #"Local Imbalance"
#set output "figures/Swap_Entanglement_LI_MI_nspin".nspin."_nspinA".nspin_A.".png"
#plot file u nspin_A:n_lim+1


set xrange [-0.1:nspin_A*log(2)+0.1]
set yrange [-0.1:1.1]
set xlabel "Mutual Information, L_A = ".nspin_A
set ylabel "IPR"
#set output "figures/Swap_Entanglement_IPR_MI_nspin".nspin."_nspinA".nspin_A.".png"
#plot file u nspin_A:n_lim+2


unset yrange
set yrange [-0.1:nspin/2*log(2)+0.1]
set xrange [-0.1:nspin_A*log(2)+0.1]
set xlabel "Mutual Information, L_A = ".nspin_A
set ylabel "Max Bipartite Entropy"
#set output "figures/Swap_Entanglement_MBE_MI_nspin".nspin."_nspinA".nspin_A.".png"
#plot file u nspin_A:n_lim+3


set yrange [-0.1:1.1]
set xrange [-0.1:nspin_A*log(2)+0.1]
set xlabel "Mutual Information, L_A = ".nspin_A
set ylabel "Imbalance Square"
#set output "figures/Swap_Entanglement_I2_MI_nspin".nspin."_nspinA".nspin_A.".png"
#plot file u nspin_A:n_lim+4

unset multiplot






########## Varying nspin/L


list_L = "6 8 10 12 14"
num_L = words(list_L)
n_lim = 2 #min(nspin/2, 2)
period="1.00"
J = "0.05"
V = "0.50"
hz = "30.00"
kick = "0.00"

#file_start="data/eigen/Sz0_DENSE_SWAP_hz_Disorder_Entanglement_nspin"
#file_end="_period".period."_iterations1_J".J."_V".V."_hz".hz."_kick".kick.".txt"
#figure_LI_IPR = "figures/Swap_Entanglement_LI_vs_IPR_Uniform_V_period".period."_iterations1_J".J."_V".V."_hz".hz."_kick".kick.".png"
#figure_MI_LI = "figures/Swap_Entanglement_MI_vs_LI_Uniform_V_nspinA".nspin_A."_period".period."_iterations1_J".J."_V".V."_hz".hz."_kick".kick.".png"

file_start="data/eigen/Sz0_DENSE_SWAP_hz_V_Disorder_Entanglement_nspin"
file_end="_period".period."_iterations1_J".J."_V".V."_hz".hz."_kick".kick.".txt"
figure_LI_IPR = "figures/Swap_Entanglement_LI_vs_IPR_Disordered_V_period".period."_iterations1_J".J."_V".V."_hz".hz."_kick".kick.".png"
figure_MI_LI = "figures/Swap_Entanglement_MI_vs_LI_Disordered_V_nspinA".nspin_A."_period".period."_iterations1_J".J."_V".V."_hz".hz."_kick".kick.".png"

list_L = "8 10 12 14"
num_L = words(list_L)
n_lim = min(nspin/2,2)
period="1.00"
J = "0.05"
V = "0.50"
hz = "16.00"
kick = "0.05"


#file_start = "data/eigen/Sz0_DENSE_SWAP_hz_V_Disorder_Neel_Overlap_nspin"
file_start="data/eigen/Sz0_DENSE_SWAP_hz_V_Disorder_Entanglement_nspin"
file_end="_period".period."_iterations1_J".J."_V".V."_hz".hz."_kick".kick.".txt"
figure_MI_LI = "figures/Swap_Entanglement_MI_vs_LI_hz_V_Disorder_nspinA".nspin_A."_period".period."_iterations1_J".J."_V".V."_hz".hz."_kick".kick.".png"

plot_title = "J = ".J.", V = ".V.", hz = ".hz.", kick = ".kick.", period = ".period

#set output figure_MI_LI
unset title
set key

set terminal png size 700,400 # tikz color # standalone
set size ratio 0.4
unset margins #-1,2,0,0

#set multiplot layout 2,2 columnsfirst title plot_title scale 1,1

set xrange [-0.1:1.1]
set yrange [-0.1:1.1]
set xlabel "Local Imbalance"
set ylabel "IPR"
set xtics 0, 0.25, 1 format "%.2f"
set ytics 0, 0.5, 1 format "%.1f"
do for [i=1:num_L] {
  set title "L = ".word(list_L,i) offset 0,-2.4
  #plot file_start.word(list_L,i).file_end u n_lim+1:n_lim+2 pt 7 notitle
}


nspin_A = 2
set xrange [-0.1:nspin_A*log(2)+0.1]
set yrange [-0.1:1.1]
set xlabel "Mutual Information, L_A =".nspin_A 
set ylabel "Local Imbalance" 
set ytics format "%.1f"
set xtics 0,nspin_A*log(2)/4, nspin_A*log(2)+0.1 format "%.1f"
do for [i=1:num_L] {
  set title "L = ".word(list_L,i) offset 0,-2.4
  #plot file_start.word(list_L,i).file_end u nspin_A:n_lim+1 pt 7 notitle
}

#unset multiplot

reset session

list_L = "8 10 12 14"
num_L = words(list_L)
n_lim = 2 #min(nspin/2,2)
period="1.00"
J = "0.05"
V = "0.50"
hz = "16.00"
kick = "0.05"

init_state = "Nayak"
file_start = "data/eigen/Sz0_DENSE_SWAP_hz_V_Disorder_".init_state."_Overlap_nspin"
file_end="_period".period."_iterations1_J".J."_V".V."_hz".hz."_kick".kick.".txt"
figure_Ovrlp_LI = "figures/Swap_Entanglement_".init_state."_Overlap_vs_LI_hz_V_Disorder_period".period."_iterations1_J".J."_V".V."_hz".hz."_kick".kick.".png"

plot_title = "J = ".J.", V = ".V.", hz = ".hz.", kick = ".kick.", period = ".period."\n{/*0.8 ".init_state." initial state}"


set xrange [-0.1:0.6]
set yrange [-0.1:1.1]
set xlabel "Overlap"
set ylabel "Local Imbalance" 
set ytics format "%.1f"
set xtics format "%.1f"
set title plot_title

set key bottom right box
set output figure_Ovrlp_LI
plot for [i=1:num_L] file_start.word(list_L,i).file_end u n_lim+2:n_lim+1 title "L = ".word(list_L,i) ps 2


set output
