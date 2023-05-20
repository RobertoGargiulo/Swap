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
########## Log Pi-Gap Swap Small J wrt J #############

file = "sort_col2_Swap_LR_spectrum_new_small_J.txt"

list_L = "4 6 8 10 12"
T0 = "1.00"
V = "3.00"
hz = "16.00"
list_J = "0.00001 0.0001 0.001 0.01 0.1"
kick = "0.001"
list_alpha = "0.50 1.00 3.00 10.00"
num_L = words(list_L)
num_J = words(list_J)
num_alpha = words(list_alpha)
f(x,y) = x+num_alpha*y ### Goes through the blocks with fixed J, where x,y change alpha, L respectively

i = 0 #alpha=0.50
j = 0 #L = 4
print i, j
print word(list_alpha,i), word(list_L,j)
set table $DataSelected
plot for [j=0:num_L-1] file every :::f(i,j)::f(i,j) u 1:2:3:4:5:6:7:12:13 w table
unset table
print "Table"
print $DataSelected

#set xrange [word(list_L,1)-0.3:word(list_L,num_L)+0.3]
#set xrange [word(list_lambda,1):word(list_L,num_L)+0.3]
set key out right box 3
unset colorbox
set size ratio 0.8
set logscale x
set xtics 1e-5, 10, 1e-1 format "10^{%T}"
set xlabel "J"
set ylabel "log10(Delta)"


i = 0 #L=4
j = 0 #alpha = 0.50
do for [i=0:num_alpha-1] {

print i

set terminal pngcairo dashed font ",13"
set output "figures/Swap_Log_Gap_Difference_J_alpha".word(list_alpha,i+1)."_L".word(list_L,num_L).".png"
set title "N_d = 2^{1-L/2}20480, T0 = ".T0.", V = ".V.", kick = ".kick.",\n hz = ".hz.", alpha = ".word(list_alpha,i+1)
plot for [j=0:num_L-1] file every :::f(i,j)::f(i,j) u 2:($12/log(10)):($13/log(10)) w errorlines title "L = ".word(list_L,j+1) lc palette frac (j+0.0)/num_L pt 4

set output "figures/Swap_Log_Gap_Comparison_J_alpha".word(list_alpha,i+1)."_L".word(list_L,num_L).".png"
plot for [j=0:num_L-1] file every :::f(i,j)::f(i,j) u 2:($14/log(10)):($15/log(10)) w errorlines title "L = ".word(list_L,j+1).", Delta_{pi}" lc palette frac (j+0.0)/num_L pt 2 ,\
     for [j=0:num_L-1] "" every :::f(i,j)::f(i,j) u 2:($16/log(10)):($17/log(10)) w errorlines title "L = ".word(list_L,j+1).", Delta_0" lc palette frac (j+0.0)/num_L pt 4 ,\
#     for [j=0:num_L-1] log10(0.01*x**(word(list_L,i+1))) lw 1.2 lc palette frac (i+0.0)/num_L dt 2 title "log10(lambda^L)"

set output "figures/Swap_Half_Spectrum_Shift_Log_Gap_Difference_J_alpha".word(list_alpha,i+1)."_L".word(list_L,num_L).".png"
plot for [j=0:num_L-1] file every :::f(i,j)::f(i,j) u 2:($18/log(10)):($19/log(10)) w errorlines title "L = ".word(list_L,j+1) lc palette frac (j+0.0)/num_L pt 4

set output "figures/Swap_Half_Spectrum_Shift_Log_Gap_Comparison_J_alpha".word(list_alpha,i+1)."_L".word(list_L,num_L).".png"
plot for [j=0:num_L-1] file every :::f(i,j)::f(i,j) u 2:($20/log(10)):($21/log(10)) w errorlines title "L = ".word(list_L,j+1).", Delta_{pi}" lc palette frac (j+0.0)/num_L pt 2 ,\
     for [j=0:num_L-1] "" every :::f(i,j)::f(i,j) u 2:($22/log(10)):($23/log(10)) w errorlines title "L = ".word(list_L,j+1).", Delta_0" lc palette frac (j+0.0)/num_L pt 4 ,\
#     for [j=0:num_L-1] log10(0.01*x**(word(list_L,i+1))) lw 1.2 lc palette frac (i+0.0)/num_L dt 2 title "log10(lambda^L)"

}

reset session

#################################################################################
########## Log Pi-Gap Swap Small J wrt L #############

file = "sort_col1_Swap_LR_spectrum_new_small_J.txt"

list_L = "4 6 8 10 12"
T0 = "1.00"
V = "3.00"
hz = "16.00"
list_J = "0.00001 0.0001 0.001 0.01 0.1"
kick = "0.001"
list_alpha = "0.50 1.00 3.00 10.00"
num_L = words(list_L)
num_J = words(list_J)
num_alpha = words(list_alpha)
f(x,y) = x+num_alpha*y ### Goes through the blocks with fixed L, where x,y change alpha, J respectively

i = 0 #alpha=0.50
j = 0 #J = 0.00001
set table $DataSelected
plot for [j=0:num_J-1] file every :::f(i,j)::f(i,j) u 1:2:3:4:5:6:7:12:13 w table
unset table
print "Table"
print $DataSelected

#set xrange [word(list_L,1)-0.3:word(list_L,num_L)+0.3]
#set xrange [word(list_lambda,1):word(list_L,num_L)+0.3]
set key out right box 3
unset colorbox
set size ratio 0.8
set xtics word(list_L,1), 2, word(list_L,num_L)
set xlabel "L"
set ylabel "log10(Delta)"


i = 0 #alpha=0.50
j = 0 #J = 0.00001
do for [i=0:num_alpha-1] {

print i

set terminal pngcairo dashed font ",13"
set output "figures/Swap_Log_Gap_Difference_L_alpha".word(list_alpha,i+1)."_L".word(list_L,num_L).".png"
set title "N_d = 2^{1-L/2}20480, T0 = ".T0.", V = ".V.", kick = ".kick.",\n hz = ".hz.", alpha = ".word(list_alpha,i+1)
plot for [j=0:num_J-1] file every :::f(i,j)::f(i,j) u 1:($12/log(10)):($13/log(10)) w errorlines title "J = ".word(list_J,j+1) lc palette frac (j+0.0)/num_J pt 4


set output "figures/figure.png"  #"figures/Swap_Log_Gap_Comparison_L_alpha".word(list_alpha,i+1)."_L".word(list_L,num_L).".png"
plot for [j=0:num_J-1] file every :::f(i,j)::f(i,j) u 1:($14/log(10)):($15/log(10)) w errorlines title "J = ".word(list_J,j+1).", Delta_{pi}" lc palette frac (j+0.0)/num_J pt 2 ,\
     for [j=0:num_J-1] "" every :::f(i,j)::f(i,j) u 1:($16/log(10)):($17/log(10)) w errorlines title "J = ".word(list_J,j+1).", Delta_0" lc palette frac (j+0.0)/num_J pt 4 ,\

set output "figures/Swap_Half_Spectrum_Shift_Log_Gap_Difference_L_alpha".word(list_alpha,i+1)."_L".word(list_L,num_L).".png"
plot for [j=0:num_J-1] file every :::f(i,j)::f(i,j) u 1:($18/log(10)):($19/log(10)) w errorlines title "J = ".word(list_J,j+1) lc palette frac (j+0.0)/num_J pt 4

set output "figures/Swap_Half_Spectrum_Shift_Log_Gap_Comparison_L_alpha".word(list_alpha,i+1)."_L".word(list_L,num_L).".png"
plot for [j=0:num_J-1] file every :::f(i,j)::f(i,j) u 1:($20/log(10)):($21/log(10)) w errorlines title "J = ".word(list_J,j+1).", Delta_{pi}" lc palette frac (j+0.0)/num_J pt 2 ,\
     for [j=0:num_J-1] "" every :::f(i,j)::f(i,j) u 1:($22/log(10)):($23/log(10)) w errorlines title "J = ".word(list_J,j+1).", Delta_0" lc palette frac (j+0.0)/num_J pt 4 ,\

}

reset session
