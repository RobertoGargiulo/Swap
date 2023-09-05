########## Gap Ratio Swap Small J wrt J #############

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

#do for [j=0:num_alpha-1] {
#  do for [k=0:num_kick-1] {
#print j, k
#print word(list_alpha,j+1), word(list_kick,k+1)
#set table $DataSelected
#plot for [i=0:num_L-1] file every :::f(i,j,k)::f(i,j,k) u 1:2:3:4:5:6:7:12:13 w table
#unset table
#print "Table"
#print $DataSelected
#  }
#}

set xrange [word(list_J,1)/1.1:word(list_J,num_J)*1.1*1.1**(num_L-1)]
set yrange [0.33:0.6]
set key out right box 3
unset colorbox
set size ratio 0.53
set logscale x
set xtics 1e-5, 10, 1e-0 format "$10^{%T}$"
set xlabel "$J$"
set ylabel "$\\langle r \\rangle$"


set margin 0
set key width 1.1 spacing 1.1
set bars 0.4

y0 = 0.386
y1 = 0.5295
x0 = word(list_J,1)
x1 = word(list_J,num_J)

fmt = ".tex"
#set terminal pngcairo dashed 
set terminal tikz color standalone 

do for [k=0:num_kick-1] {

kick = word(list_kick,k+1)
set output "figures/Swap_Gap_Ratio_wrt_J_kick".kick."_L".word(list_L,num_L)."".fmt
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
  set arrow from x0,y0 to x1,y0 nohead dt 2
  set arrow from x0,y1 to x1,y1 nohead dt 2
  set label "Poisson" at x1/15,y0-0.02 font ",9"
  set label "COE" at x1/15,y1-0.02 font ",9"
  plot for [i=1:num_L-1] file every :::f(i,j,k)::f(i,j,k) u ($2*1.1**(i)):8:9 w errorlines title "$L = ".word(list_L,i+1)."$" lc palette frac (i-1+0.0)/(num_L-1) pt 4 ps 0.5

  j = 1
  print j, k
  alpha = word(list_alpha,j+1) 
  print alpha, kick
  set key at screen 0.6, 0. width 1 spacing 1.1 font ",9" maxrows 1 #inside top left 
  set xtics 1e-5, 10, 1e-0 format "$10^{%T}$"
  set xlabel "$J$"
  set title "(b) $\\alpha = ".alpha."$"
  plot for [i=1:num_L-1] file every :::f(i,j,k)::f(i,j,k) u ($2*1.1**(i)):8:9 w errorlines title "$L = ".word(list_L,i+1)."$" lc palette frac (i-1+0.0)/(num_L-1) pt 4 ps 0.5
}
unset multiplot

reset session
