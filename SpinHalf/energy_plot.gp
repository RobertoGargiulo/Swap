set terminal png # tikz color # standalone
#set output "figures/figure.png"
set size ratio 0.5
set margins -1,2,0,0


file="data/eigen/Sz0_DENSE_SWAP_hz_Disorder_Quasi_Energie_nspin10_period1.00_iterations1_J0.05_V0.50_hz16.00_kick0.00.txt"

dim=system(sprintf("wc -l %s | awk '{print $1}'",file))
#print dim

L = "10"
J = "0.05"
V = "0.50"
hz = "16.00"
kick = "0.00"
period = "1.00"

set ytics -pi, pi/2, pi format "%.2f"
set xrange [0:dim/2-1]
set yrange [-pi-0.1:pi+0.1]

set title "Single Iteration, nspin = ".L.", J = ".J.", V = ".V.", hz = ".hz.", kick = ".kick.", period = ".period
set output "figures/Energy_Gap.png"
plot file every ::0::dim/2-1 notitle, "" every ::dim/2::dim-1 notitle
