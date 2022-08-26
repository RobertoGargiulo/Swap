set terminal png # tikz color # standalone
set size ratio 0.5
set margins -1,2,0,0

max(x, y) = (x > y ? x : y)
min(x, y) = (x < y ? x : y)


###### Parameters
nspin = "12"
n_lim = min(nspin/2, 4)
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




set yrange [-0.1:2.1]
set xrange [-0.1:nspin_A*log(2)+0.1]
set xlabel "Mutual Information L_A =".nspin_A 
set ylabel "Local Imbalance" 
set output "figures/figure.png"
plot file1 u nspin_A:n_lim+1 title "J = 0.00", file2 u nspin_A:n_lim+1 title "J = 0.05"



plot_title = "nspin = ".nspin.", J = ".J.", V = ".V.", hz = ".hz.", kick = ".kick.", period = ".period

#set output "figures/Swap_Entanglement_MI_corr".figure_param
unset title
set xtics 0, nspin_A*log(2)/3, nspin_A*log(2) format "%.1f"

set multiplot layout 2,2 columnsfirst scale 1,1 title plot_title

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

set output
