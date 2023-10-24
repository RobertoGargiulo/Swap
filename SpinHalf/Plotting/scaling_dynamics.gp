reset session
########## sigmaz Dynamics wrt L #############

init_state = "Neel"
list_L = "4 6 8 10 12"
list_dis = "10240 5120 2560 1280 640"
T0 = "1.00"
V = "3.00"
hz = "16.00"
steps = "10000"
list_J = "0.00010 0.00100 0.01000 0.10000 1.00000"
list_kick = "0.000 0.010"
list_alpha = "0.50 3.00"
num_L = words(list_L)
num_J = words(list_J)
num_kick = words(list_kick)
num_alpha = words(list_alpha)

#do for [i=0:num_J-1] {
#  do for [j=0:num_alpha-1] {
#    do for [k=0:num_kick-1] {
#  
#print i, j, k
#Jxy = word(list_J,i+1)
#alpha = word(list_alpha,j+1)
#kick = word(list_kick,k+1)
#print Jxy, " ", alpha, " ", kick
#
#print "Table"
#
#set table $DataSelected
#
#plot for [q=0:num_L-1] "data/dynamics/sigmaz_Swap_LR_".init_state."_nspin".word(list_L,q+1)."_period".T0."_n_disorder".word(list_dis,q+1)."_Jxy".Jxy."_Vzz".V."_hz".hz."_kick".kick."_alpha".alpha.".txt" \
#every ::1::17 u 1:( column(1 + word(list_L,q+1)/4 * 2) ) : ( column(1 + word(list_L,q+1) - 1) ) w table
#unset table
#print $DataSelected
#
#    }
#  }
#}

#set size ratio 0.5
#set margins -1,2,0,0
set key out right box 3
unset colorbox
#set size ratio 0.8
set logscale x
set yrange [-1.1:1.1]
set xlabel "t"
set ylabel "\\sigma^z_{2(L/4)}(t)"

do for [i=0:num_J-1] {
  do for [j=0:num_alpha-1] {
    do for [k=0:num_kick-1] {
  
print i, j, k
Jxy = word(list_J,i+1)
alpha = word(list_alpha,j+1)
kick = word(list_kick,k+1)
print Jxy, " ", alpha, " ", kick
do for [q=0:num_L-1] {
  print "columns: ", 1, 1 + word(list_L,q+1)/4 * 2, word(list_L,q+1)
}

#set terminal pngcairo dashed font ",13"
#set output "figures/Swap_Dynamics_".init_state."_sigmaz_avg_wrt_L_J".Jxy."_kick".kick."_alpha".alpha."_L".word(list_L,num_L).".png"
#set title "N_d = 2^{1-L/2}20480, steps = ".steps.", T0 = ".T0.",\n J = ".Jxy.", V = ".V.", hz = ".hz.",\n kick = ".kick.", alpha = ".alpha
#plot for [q=0:num_L-1] "data/dynamics/sigmaz_Swap_LR_".init_state."_nspin".word(list_L,q+1)."_period".T0."_n_disorder".word(list_dis,q+1)."_Jxy".Jxy."_Vzz".V."_hz".hz."_kick".kick."_alpha".alpha.".txt" \
#  u 1:( (-1)**($1/T0) * column( 1 + word(list_L,q+1)/4 * 2  ) ) w lp title "L = ".word(list_L,q+1) lc palette frac (q+0.0)/num_L


#set output "figures/Swap_Dynamics_sigmaz_err_wrt_L_J".Jxy."kick".kick."_alpha".alpha."_L".word(list_L,num_L).".png"
#set title "N_d = 2^{1-L/2}20480, T0 = ".T0.", V = ".V.", kick = ".kick.",\n hz = ".hz.", alpha = ".alpha
#plot for [i=0:num_L-1] file u 1:column( 1 + word(list_L,q+1)/4 * 2  ) w errorlines title "L = ".word(list_L,i+1) lc palette frac (i+0.0)/num_L pt 4
  
    }
  }
}

#set terminal png
#set output figure
#plot for [i=1:4] "data/dynamics/sigmaz_Swap_vs_MBL_Neel_nspin".word(list_L,i)."_steps".steps."_period".T0."_n_disorder".word(list_dis,i)."_J".Jxy."_V".V."_hz".hz."_kick".kick.".txt" \
#       u 1:(                        -(sum[k=1:int(word(list_L,i))] (-1)**k * column(k+1)) / word(list_L,i) ) w lp lw 1.6 title "MBL, L = ".word(list_L,i) ,\
#    for [i=1:4] "" u 1:(  (-1)**($1/T0) * (sum[k=1:int(word(list_L,i))] (-1)**k * column(k+word(list_L,i)+1)) / word(list_L,i) ) w lp lw 1.6 title "Swap, L = ".word(list_L,i)  


reset session




init_state = "HalfNeel"
list_L = "4 6 8 10 12"
list_dis = "10240 5120 2560 1280 640"
T0 = "1.00"
V = "3.00"
hz = "16.00"
steps = "10000"
list_J = "0.00010 0.00100 0.01000 0.10000 1.00000"
list_kick = "0.000 0.010"
list_alpha = "0.50 3.00"
num_L = words(list_L)
num_J = words(list_J)
num_kick = words(list_kick)
num_alpha = words(list_alpha)


#do for [i=0:num_J-1] {
#  do for [j=0:num_alpha-1] {
#    do for [k=0:num_kick-1] {
#  
#print i, j, k
#Jxy = word(list_J,i+1)
#alpha = word(list_alpha,j+1)
#kick = word(list_kick,k+1)
#print Jxy, " ", alpha, " ", kick
#
#print "Table"
#
#
#do for [q=0:num_L-1] {
#  set table $DataSelected
#  filename= "data/dynamics/sigmaz_Swap_LR_".init_state."_nspin".word(list_L,q+1)."_period".T0."_n_disorder".word(list_dis,q+1)."_Jxy".Jxy."_Vzz".V."_hz".hz."_kick".kick."_alpha".alpha.".txt"
#  sigmaz_0=system("head -9 ".filename." | tail -1")
#  print "L = ", word(list_L,q+1), "sigmaz_0 = ", sigmaz_0
#  plot filename \
#  every ::10000::10020 u 1:2:3:((sum[k=1:int(word(list_L,q+1)/2)] sgn( int(word(sigmaz_0,1+2*k-1) - word(sigmaz_0,1+2*k))  )  * (column(1+2*k-1) - column(1+2*k))) / int(word(list_L,q+1))) w table
#  unset table
#  print $DataSelected
#}
#
#    }
#  }
#}


set terminal pngcairo dashed font ",13"


array time[3,steps]
do for [i=1:steps] { 
  time[1,i] = i
  time[2,i] = 2*i
  time[3,i] = 3*i
}
set output "figures/figure.png"
plot sample [i=1:steps] '+' u (time[1,i]):(time[2,i]) w lp pt 1,\
     sample [i=1:steps] '+' u (time[1,i]):(time[3,i]) w lp pt 2

set key out right box 3
unset colorbox
#set size ratio 0.8
#set logscale x
#set yrange [-1.1:1.1]
set xlabel "t"
set ylabel "\\sigma^z_{2(L/4)}(t)"

set terminal pngcairo dashed font ",13"

do for [i=0:num_J-1] {
  do for [j=0:num_alpha-1] {
    do for [k=0:num_kick-1] {
  
print i, j, k
Jxy = word(list_J,i+1)
alpha = word(list_alpha,j+1)
kick = word(list_kick,k+1)
print Jxy, " ", alpha, " ", kick


set output "figures/Swap_Dynamics_".init_state."_Z_avg_wrt_L_J".Jxy."_kick".kick."_alpha".alpha."_L".word(list_L,num_L).".png"
do for [q=0:num_L-1] {
  filename= "data/dynamics/sigmaz_Swap_LR_".init_state."_nspin".word(list_L,q+1)."_period".T0."_n_disorder".word(list_dis,q+1)."_Jxy".Jxy."_Vzz".V."_hz".hz."_kick".kick."_alpha".alpha.".txt"
  sigmaz_0=system("head -9 ".filename." | tail -1 | awk '{for(i=2;i<=1+".word(list_L,q+1).";i++) printf $i\" \"; print \"\"}'")
  imb_0=system("echo ".sigmaz_0." | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {NUM=0; for(i=1;i<=NF/2;i++){NUM=NUM+abs($(2*i-1)-$(2*i))}; printf \"%.17g\\n\", NUM/".word(list_L,q+1)."}' ")
  imb_0 = imb_0 + 0.0
  print "L = ", word(list_L,q+1), "\n sigmaz_0 = ", sigmaz_0, "\n imb_0 = ", imb_0, "\n"
  plot filename \
    u 1:((-1)**($1-1)*(sum[k=1:int(word(list_L,q+1)/2)] sgn( (word(sigmaz_0,2*k-1)+0.0) - (word(sigmaz_0,2*k)+0.0)  )  * (column(1+2*k-1) - column(1+2*k)))/ int(word(list_L,q+1))/imb_0) \
    w lp title "L = ".word(list_L,q+1) lc palette frac (q+0.0)/num_L
}

    }
  }
}




reset session


########################################

