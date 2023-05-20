########## Log Pi-Gap Khemani wrt lambda #############

file = "Khemani_spectrum_L10_kick0.txt"

# In generale: i + num_J * j + num_J * num_V * k, 0<i<num_J-1; 0<j<num_V-1; 0<k<num_hz-1

list_L = "4 6 8 10"
list_Nd = "10240 5120 2560 1280"
T0 = "1.00"
V = "1.00"
list_lambda = "0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28 2.56"
kick = "0.00"
num_L = words(list_L)
num_lambda = words(list_lambda)

#set xrange [word(list_L,1)-0.3:word(list_L,num_L)+0.3]
#set xrange [word(list_lambda,1):word(list_L,num_L)+0.3]
set key out right box 3
unset colorbox
set size ratio 0.8
set logscale x
set xtics 0.01, 4, 2.56
set xlabel "$lambda$"
set ylabel "log10(Delta)"

set terminal pngcairo dashed font ",13"
set output "figures/Khemani_Log_Gap_Difference_lambda_L".word(list_L,num_L).".png"
set title "N_d = 2^{1-L/2}20480, T0 = ".T0.", V = ".V.", kick = ".kick.",\n hx = lambda * 0.1, hy = lambda * 0.15, hz = lambda * 0.45"
plot for [i=0:num_L-1] file every :::i::i u 2:($13/log(10)):($14/log(10)) w errorlines title "L = ".word(list_L,i+1) lc palette frac (i+0.0)/num_L pt 4

set output "figures/Khemani_Log_Gap_Comparison_lambda_L".word(list_L,num_L).".png"
plot for [i=0:num_L-1] file every :::i::i u 2:($15/log(10)):($16/log(10)) w errorlines title "L = ".word(list_L,i+1).", Delta_{pi}" lc palette frac (i+0.0)/num_L pt 2 ,\
     for [i=0:num_L-1] "" every :::i::i u 2:($17/log(10)):($18/log(10)) w errorlines title "L = ".word(list_L,i+1).", Delta_0" lc palette frac (i+0.0)/num_L pt 4 ,\
#     for [i=0:num_L-1] log10(0.01*x**(word(list_L,i+1))) lw 1.2 lc palette frac (i+0.0)/num_L dt 2 title "log10(lambda^L)"

set output "figures/Khemani_Half_Spectrum_Shift_Log_Gap_Difference_lambda_L".word(list_L,num_L).".png"
set title "N_d = 2^{1-L/2}20480, T0 = ".T0.", V = ".V.", kick = ".kick.",\n hx = lambda * 0.1, hy = lambda * 0.15, hz = lambda * 0.45"
plot for [i=0:num_L-1] file every :::i::i u 2:($19/log(10)):($20/log(10)) w errorlines title "L = ".word(list_L,i+1) lc palette frac (i+0.0)/num_L pt 4

set output "figures/Khemani_Half_Spectrum_Shift_Log_Gap_Comparison_lambda_L".word(list_L,num_L).".png"
plot for [i=0:num_L-1] file every :::i::i u 2:($21/log(10)):($22/log(10)) w errorlines title "L = ".word(list_L,i+1).", Delta_{pi}" lc palette frac (i+0.0)/num_L pt 2 ,\
     for [i=0:num_L-1] "" every :::i::i u 2:($23/log(10)):($24/log(10)) w errorlines title "L = ".word(list_L,i+1).", Delta_0" lc palette frac (i+0.0)/num_L pt 4 ,\

reset session

#################################################################################
########## Log Pi-Gap Khemani wrt L #############

file = "sort_col1_Khemani_spectrum_L10_kick0.txt"

# In generale: i + num_J * j + num_J * num_V * k, 0<i<num_J-1; 0<j<num_V-1; 0<k<num_hz-1

list_L = "4 6 8 10"
list_Nd = "10240 5120 2560 1280"
T0 = "1.00"
V = "1.00"
list_lambda = "0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28 2.56"
kick = "0.00"
num_L = words(list_L)
num_lambda = words(list_lambda)

#set xrange [word(list_L,1)-0.3:word(list_L,num_L)+0.3]
#set xrange [word(list_lambda,1):word(list_L,num_L)+0.3]
#set key out right box 3
unset key
unset colorbox
set size ratio 0.8
set xtics word(list_L,1), 2, word(list_L,num_L)
set xlabel "L"
set ylabel "log10(Delta)"

set terminal pngcairo dashed font ",13"
set output "figures/Khemani_Log_Gap_Difference_wrt_L_L".word(list_L,num_L).".png"
set title "N_d = 2^{1-L/2}20480, T0 = ".T0.", V = ".V.", kick = ".kick.",\n hx = lambda * 0.1, hy = lambda * 0.15, hz = lambda * 0.45"
plot for [i=0:num_lambda-1] file every :::i::i u 1:($13/log(10)):($14/log(10)) w errorlines lc palette frac (i+0.0)/num_lambda pt 4

set output "figures/Khemani_Log_Gap_Comparison_wrt_L_L".word(list_L,num_L).".png"
plot for [i=0:num_lambda-1] file every :::i::i u 1:($15/log(10)):($16/log(10)) w errorlines lc palette frac (i+0.0)/num_lambda pt 2 ,\
     for [i=0:num_lambda-1] "" every :::i::i u 1:($17/log(10)):($18/log(10)) w errorlines lc palette frac (i+0.0)/num_lambda pt 4 ,\

set output "figures/Khemani_Half_Spectrum_Shift_Log_Gap_Difference_wrt_L_L".word(list_L,num_L).".png"
set title "N_d = 2^{1-L/2}20480, T0 = ".T0.", V = ".V.", kick = ".kick.",\n hx = lambda * 0.1, hy = lambda * 0.15, hz = lambda * 0.45"
plot for [i=0:num_lambda-1] file every :::i::i u 1:($19/log(10)):($20/log(10)) w errorlines lc palette frac (i+0.0)/num_lambda pt 4

set output "figures/Khemani_Half_Spectrum_Shift_Log_Gap_Comparison_wrt_L_L".word(list_L,num_L).".png"
plot for [i=0:num_lambda-1] file every :::i::i u 1:($21/log(10)):($22/log(10)) w errorlines lc palette frac (i+0.0)/num_lambda pt 2 ,\
     for [i=0:num_lambda-1] "" every :::i::i u 1:($23/log(10)):($24/log(10)) w errorlines lc palette frac (i+0.0)/num_lambda pt 4 ,\

reset session

#################################################################################
