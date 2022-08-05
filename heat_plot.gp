set terminal png # tikz color # standalone
#set output "figures/figure.png"
set size ratio 1
#set margins -1,1,0,0

set cbrange [0.30:0.57]

set ytics format "%.1f"
set xtics format "%.1f"
set ylabel "h_{z}"
set xlabel "V"
set cblabel "Average Gap Ratio"


nspin="12"
iter="80"
J="0.50"
set title "no kick, Sz=0 subspace, n_{spin} = ".nspin.", n_{iter} = ".iter.", J = ".J
filename="data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_Gap_Ratio_nspin".nspin."_iterations".iter."_J".J."_up_to_V3_hz6.txt"

#set xrange [0:3]
#set yrange [0:6]
#set output "figures/MBL_Full_PD_Gap_Ratio_nspin".nspin."_J".J."_Dense_up_toV3_hz6.png"

set xrange [0:0.75]
set yrange [0:3]
set output "figures/MBL_Full_PD_Gap_Ratio_nspin".nspin."_J".J."_Dense_up_toV0.75_hz3.png"


set view map
#set dgrid3d 100,100,2
set pm3d at b
unset key
unset surface
splot filename u 2:3:4

set output
