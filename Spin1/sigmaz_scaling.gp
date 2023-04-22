reset session

set size ratio 0.5
set margins -1,2,0,0

n_dis2 = 1280

L = "8"
steps = "10000"
T0 = "1.00"
n_dis = "160"
J = "0.60"
V = "2.70"
hz = "4.00"
kick = "0.00"

Data = "data/dynamics/sigmaz_MBL_Neel_nspin".L."_steps".steps."_period".T0."_n_disorder".n_dis."_J".J."_V".V."_hz".hz.".txt"
figure = "figures/sigmaz_MBL_Neel_nspin".L."_steps".steps."_period".T0."_n_disorder".n_dis."_J".J."_V".V."_hz".hz.".png"
f_title = "L = ".L.", steps = ".steps.", T0 = ".T0.", n_{dis} = ".n_dis.",\n J = ".J.", V = ".V.", hz = ".hz

set title f_title

set logscale x
set yrange [0:1]
set xlabel "t"
set ylabel "I(t)"

#set terminal png
#set output figure
#plot Data u 1:( -(sum[i=1:int(L)] (-1)**i * column(i+1))/L ) w lp lw 1.6 title "MBL" ,\


########################################################################

reset session

set size ratio 0.5
set margins -1,2,0,0

L = "8"
steps = "10000"
T0 = "1.00"
n_dis = "320"
J = "0.15"
V = "3.50"
hz = "8.00"
kick = "0.10"

Data = "data/dynamics/sigmaz_Swap_vs_MBL_Neel_nspin".L."_steps".steps."_period".T0."_n_disorder".n_dis."_J".J."_V".V."_hz".hz."_kick".kick.".txt"
figure = "figures/sigmaz_Swap_vs_MBL_Neel_nspin".L."_steps".steps."_period".T0."_n_disorder".n_dis."_J".J."_V".V."_hz".hz."_kick".kick.".png"
f_title = "L = ".L.", steps = ".steps.", T0 = ".T0.", n_{dis} = ".n_dis.",\n J = ".J.", V = ".V.", hz = ".hz.", kick = ".kick

set title f_title

set logscale x
set yrange [0:1]
set xlabel "t"
set ylabel "I(t)"

#set terminal png
#set output figure
#plot Data u 1:( -(sum[i=1:int(L)] (-1)**i * column(i+1))/L ) w lp lw 1.6 title "MBL" ,\
#     Data u 1:( (-1)**($1/T0) * (sum[i=1:int(L)] (-1)**i * column(i+L+1))/L ) w lp lw 1.6 title "Swap" ,\


########################################################################

reset session

set size ratio 0.5
set margins -1,2,0,0

steps = "10000"
T0 = "1.00"
J = "0.15"
V = "3.50"
hz = "8.00"
kick = "0.10"

figure = "figures/sigmaz_MBL_Neel_compare_nspin_steps".steps."_period".T0."_J".J."_V".V."_hz".hz."_kick".kick.".png"
f_title = "steps = ".steps.", T0 = ".T0.", n_{dis} = 2560*2^{1-L/2},\n J = ".J.", V = ".V.", hz = ".hz

set title f_title

set logscale x
#set logscale y
set yrange [0.8:1.1]
set xlabel "t"
set ylabel "I(t)"

list_L = "2 4 6 8"
list_dis = "2560 1280 640 320"

#set terminal png
#set output figure
#plot for [i=1:4] "data/dynamics/sigmaz_Swap_vs_MBL_Neel_nspin".word(list_L,i)."_steps".steps."_period".T0."_n_disorder".word(list_dis,i)."_J".J."_V".V."_hz".hz."_kick".kick.".txt" \
#       u 1:(                -(sum[k=1:int(word(list_L,i))] (-1)**k * column(k+1))                    /word(list_L,i) ) w lp lw 1.6 title "MBL, L = ".word(list_L,i) ,\
#    "" u 1:( (-1)**($1/T0) * (sum[k=1:2] (-1)**k * column(k+2+1)) / word(list_L,i) ) w lp lw 1.6 title "Swap, L = ".word(list_L,i) ,\
#int(word(list_L,i))




########################################################################

reset session

set size ratio 0.5
set margins -1,2,0,0

steps = "10000"
T0 = "1.00"
J = "0.15"
V = "3.50"
hz = "8.00"
kick = "0.10"

figure = "figures/sigmaz_Swap_Neel_compare_nspin_steps".steps."_period".T0."_J".J."_V".V."_hz".hz."_kick".kick.".png"
f_title = "steps = ".steps.", T0 = ".T0.", n_{dis} = 2560*2^{1-L/2},\n J = ".J.", V = ".V.", hz = ".hz.", kick = ".kick

set title f_title

set logscale x
#set logscale y
set yrange [-1.1:1.1]
set xlabel "t"
set ylabel "I(t)"

list_L = "2 4 6 8"
list_dis = "2560 1280 640 320"

#set terminal png
#set output figure
#plot for [i=1:4] "data/dynamics/sigmaz_Swap_vs_MBL_Neel_nspin".word(list_L,i)."_steps".steps."_period".T0."_n_disorder".word(list_dis,i)."_J".J."_V".V."_hz".hz."_kick".kick.".txt" \
#     u 1:( (-1)**(column(1)/T0) * (sum[k=1:int(word(list_L,i))] (-1)**k * column(int(k+word(list_L,i)+1))) / word(list_L,i) ) w lp lw 1.6 title "Swap, L = ".word(list_L,i) ,\


########################################################################

reset session

set size ratio 0.5
set margins -1,2,0,0

steps = "10000"
T0 = "1.00"
J = "0.15"
V = "3.50"
hz = "8.00"
kick = "0.10"

figure = "figures/sigmaz_Swap_vs_MBL_Neel_compare_nspin_steps".steps."_period".T0."_J".J."_V".V."_hz".hz."_kick".kick.".png"
f_title = "steps = ".steps.", T0 = ".T0.", n_{dis} = 2560*2^{1-L/2},\n J = ".J.", V = ".V.", hz = ".hz.", eps = ".kick

set title f_title

set logscale x
#set logscale y
set yrange [0.8:1.1]
set xlabel "t"
set ylabel "I(t)"

list_L = "2 4 6 8"
list_dis = "2560 1280 640 320"

set terminal png
set output figure
plot for [i=1:4] "data/dynamics/sigmaz_Swap_vs_MBL_Neel_nspin".word(list_L,i)."_steps".steps."_period".T0."_n_disorder".word(list_dis,i)."_J".J."_V".V."_hz".hz."_kick".kick.".txt" \
       u 1:(                        -(sum[k=1:int(word(list_L,i))] (-1)**k * column(k+1)) / word(list_L,i) ) w lp lw 1.6 title "MBL, L = ".word(list_L,i) ,\
    for [i=1:4] "" u 1:(  (-1)**($1/T0) * (sum[k=1:int(word(list_L,i))] (-1)**k * column(k+word(list_L,i)+1)) / word(list_L,i) ) w lp lw 1.6 title "Swap, L = ".word(list_L,i)  

#








