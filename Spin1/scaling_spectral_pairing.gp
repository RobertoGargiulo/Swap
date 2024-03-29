########## Log Pi-Gap Difference (Compact Multiplot) Swap wrt J #############

file = "sort_col2_Swap_LR_spectrum.txt"

list_L = "2 4 6 8"
T0 = "1.00"
V = "3.00"
hz = "16.00"
list_J = "0.00001 0.0001 0.001 0.01 0.1 1.0"
list_kick = "0.0"
list_alpha = "0.50 3.00"
num_L = words(list_L)
num_J = words(list_J)
num_kick = words(list_kick)
num_alpha = words(list_alpha)
f(x,y,z) = x+num_L*y+num_alpha*num_L*z ### Goes through the blocks with varying J, where x,y,z change alpha, L, kick respectively

set xrange [word(list_J,1):word(list_J,num_J)]
set key out right box 3
unset colorbox
set size ratio 0.53
set logscale x
set xtics 1e-5, 10, 1e-0 format "$10^{%T}$"
set xlabel "$J$"
set ylabel "$\\ell$"


set margin 0
set key width 1.1 spacing 1.1
set bars 0.4

fmt = ".tex"
#set terminal pngcairo dashed 
set terminal tikz color standalone 

do for [k=0:num_kick-1] {

kick = word(list_kick,k+1)
set output "figures/Swap_Log_Gap_Difference_wrt_J_kick".kick."_L".word(list_L,num_L)."".fmt
set multiplot layout num_alpha,1 columns margins 0.12,1,0.12,0.95 \
  spacing 0,0.09 

  j = 0
  print j, k
  alpha = word(list_alpha,j+1) 
  print alpha, kick
  unset key
  set xtics format ""
  unset xlabel
  set title "(a) $\\alpha = ".alpha."$"
  plot for [i=0:num_L-1] file every :::f(i,j,k)::f(i,j,k) u 2:($12/log(10)):($13/log(10)) w errorlines title "$L = ".word(list_L,i+1)."$" lc palette frac (i+0.0)/num_L pt 4

  j = 1
  print j, k
  alpha = word(list_alpha,j+1) 
  print alpha, kick
  set key width 1 spacing 1.1 font ",8" inside bottom right
  set xtics 1e-5, 10, 1e-0 format "$10^{%T}$"
  set xlabel "$J$"
  set title "(b) $\\alpha = ".alpha."$"
  plot for [i=0:num_L-1] file every :::f(i,j,k)::f(i,j,k) u 2:($12/log(10)):($13/log(10)) w errorlines title "$L = ".word(list_L,i+1)."$" lc palette frac (i+0.0)/num_L pt 4
}
unset multiplot

reset session

#################################################################################

########## Log Pi-Gap Swap Small J wrt J #############

file = "sort_col2_Swap_LR_spectrum.txt"

list_L = "2 4 6 8"
T0 = "1.00"
V = "3.00"
hz = "16.00"
list_J = "0.00001 0.0001 0.001 0.01 0.1 0.5"
list_kick = "0.0"
list_alpha = "0.50 3.00"
num_L = words(list_L)
num_J = words(list_J)
num_kick = words(list_kick)
num_alpha = words(list_alpha)
f(x,y,z) = x+num_L*y+num_alpha*num_L*z ### Goes through the blocks with varying J, where x,y,z change alpha, L, kick respectively

#do for [j=0:num_alpha-1] {
#  do for [k=0:num_kick-1] {
##print j, k
#print "alpha = ", word(list_alpha,j+1), ", kick = ", word(list_kick,k+1)
#set table $DataSelected
#plot for [i=0:num_L-1] file every :::f(i,j,k)::f(i,j,k) u 1:2:3:4:5:6:7:12:13 w table
#unset table
#print "Table"
#print $DataSelected
#  }
#}

set key out right box 3
unset colorbox
set size ratio 0.8
set logscale x
set xtics 1e-5, 10, 1e-0 format "$10^{%T}$"
set xlabel "$J$"
set ylabel "$\\log_{10}(\\Delta)$"


#do for [j=0:num_alpha-1] {
#  do for [k=0:num_kick-1] {
#
#print j, k
#alpha = word(list_alpha,j+1) 
#kick = word(list_kick,k+1)
#print alpha, kick
#
#set margin 0
#set key width 1.1 spacing 1.1
#
#
#fmt = ".tex"
##set terminal pngcairo dashed font ",13"
#set terminal tikz color standalone scale 0.65, 0.65
#
#set output "figures/Swap_Log_Gap_Difference_wrt_J_kick".kick."_alpha".alpha."_L".word(list_L,num_L)."".fmt
##set title "N_d = 2^{1-L/2}20480, T0 = ".T0.", V = ".V.", kick = ".kick.",\n hz = ".hz.", alpha = ".alpha
#set notitle
#plot for [i=0:num_L-1] file every :::f(i,j,k)::f(i,j,k) u 2:($12/log(10)):($13/log(10)) w errorlines title "$L = ".word(list_L,i+1)."$" lc palette frac (i+0.0)/num_L pt 4
#
#set output "figures/Swap_Log_Gap_Comparison_wrt_J_kick".kick."_alpha".alpha."_L".word(list_L,num_L)."".fmt
#plot for [i=0:num_L-1] file every :::f(i,j,k)::f(i,j,k) u 2:($14/log(10)):($15/log(10)) w errorlines  lc palette frac (i+0.0)/num_L pt 2 notitle , \
#     for [i=0:num_L-1] "" every :::f(i,j,k)::f(i,j,k) u 2:($16/log(10)):($17/log(10)) w errorlines lc palette frac (i+0.0)/num_L pt 4 notitle, \
#     for [i=0:num_L-1] NaN title "$L = ".word(list_L,i+1)."$" lc palette frac (i+0.0)/num_L, NaN w p lc "black" pt 2 title "$\\Delta_{\\pi}$", NaN w p lc "black" pt 4 title "$\\Delta_0$"
#     
#
#set output "figures/Swap_Half_Spectrum_Shift_Log_Gap_Difference_wrt_J_kick".kick."_alpha".alpha."_L".word(list_L,num_L)."".fmt
#plot for [i=0:num_L-1] file every :::f(i,j,k)::f(i,j,k) u 2:($18/log(10)):($19/log(10)) w errorlines title "$L = ".word(list_L,i+1)."$" lc palette frac (i+0.0)/num_L pt 4
#
#set output "figures/Swap_Half_Spectrum_Shift_Log_Gap_Comparison_wrt_J_kick".kick."_alpha".alpha."_L".word(list_L,num_L)."".fmt
#plot for [i=0:num_L-1] file every :::f(i,j,k)::f(i,j,k) u 2:($20/log(10)):($21/log(10)) w errorlines title "$L = ".word(list_L,i+1)."$".", $\\Delta_{\\pi}$" lc palette frac (i+0.0)/num_L pt 2 ,\
#     for [i=0:num_L-1] "" every :::f(i,j,k)::f(i,j,k) u 2:($22/log(10)):($23/log(10)) w errorlines title "$L = ".word(list_L,i+1)."$".", $\\Delta_0$" lc palette frac (i+0.0)/num_L pt 4
#
#  }
#}

reset session
