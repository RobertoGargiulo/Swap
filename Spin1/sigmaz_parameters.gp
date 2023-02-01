reset session #Swap/MBL varying kick strength

set size ratio 0.5
set margins -1,2,0,0

steps = "10000"
T0 = "1.00"
J = "0.05"
V = "10.00"
hz = "20.00"
L = "8"
n_dis = "320"

figure = "figures/sigmaz_Swap_vs_MBL_Neel_compare_kick_nspin".L."_steps".steps."_period".T0."_J".J."_V".V."_hz".hz.".png"
f_title = "L = ".L.", steps = ".steps.", T0 = ".T0.", n_{dis} = 2560*2^{1-L/2},\n J = ".J.", V = ".V.", hz = ".hz

set title f_title

set logscale x
#set logscale y
set yrange [-.1:1.1]
set xlabel "t"
set ylabel "I(t)"

list_kick = "0.00 0.05 0.10"
num_kick = words(list_kick)

set terminal png
set output figure
plot "data/dynamics/sigmaz_Swap_vs_MBL_Neel_nspin".L."_steps".steps."_period".T0."_n_disorder".n_dis."_J".J."_V".V."_hz".hz."_kick".word(list_kick,1).".txt" \
      u 1:( -(sum[k=1:int(L)] (-1)**k * column(k+1)) / L ) w lp lw 1.6 title "MBL" ,\
      for [i=1:num_kick] "data/dynamics/sigmaz_Swap_vs_MBL_Neel_nspin".L."_steps".steps."_period".T0."_n_disorder".n_dis."_J".J."_V".V."_hz".hz."_kick".word(list_kick,i).".txt" \
      u 1:( (-1)**(column(1)/T0) * ( sum[k=1:int(L)] (-1)**k * column(int(k+L+1)) ) / L )  w lp lw 1.6 title "Swap, kick = ".word(list_kick,i) ,\



########################################################################



reset session #Swap/MBL varying hz 

set size ratio 0.5
set margins -1,2,0,0

steps = "10000"
T0 = "1.00"
J = "0.05"
V = "10.00"
hz = "20.00"
L = "8"
n_dis = "320"

figure = "figures/sigmaz_Swap_vs_MBL_Neel_compare_kick_nspin".L."_steps".steps."_period".T0."_J".J."_V".V."_hz".hz.".png"
f_title = "L = ".L.", steps = ".steps.", T0 = ".T0.", n_{dis} = 2560*2^{1-L/2},\n J = ".J.", V = ".V.", hz = ".hz

set title f_title

set logscale x
#set logscale y
set yrange [-.1:1.1]
set xlabel "t"
set ylabel "I(t)"

list_kick = "0.00 0.05 0.10"
num_kick = words(list_kick)

set terminal png
set output figure
plot "data/dynamics/sigmaz_Swap_vs_MBL_Neel_nspin".L."_steps".steps."_period".T0."_n_disorder".n_dis."_J".J."_V".V."_hz".hz."_kick".word(list_kick,1).".txt" \
      u 1:( -(sum[k=1:int(L)] (-1)**k * column(k+1)) / L ) w lp lw 1.6 title "MBL" ,\
      for [i=1:num_kick] "data/dynamics/sigmaz_Swap_vs_MBL_Neel_nspin".L."_steps".steps."_period".T0."_n_disorder".n_dis."_J".J."_V".V."_hz".hz."_kick".word(list_kick,i).".txt" \
      u 1:( (-1)**(column(1)/T0) * ( sum[k=1:int(L)] (-1)**k * column(int(k+L+1)) ) / L )  w lp lw 1.6 title "Swap, kick = ".word(list_kick,i) ,\



