########## Log Pi-Gap #############

file = "sort_col1_Swap_LR_spectrum_L12_V_hz_alpha_large_N_dis.txt"

# In generale: i + num_J * j + num_J * num_V * k, 0<i<num_J-1; 0<j<num_V-1; 0<k<num_hz-1

list_L = "4 6 8 10"
list_Nd = "10240 5120 2560 1280"
T0 = "1.00"
J = "0.10"
list_V = "1.00 2.00 3.00"
list_hz = "2.00 4.00 8.00 16.00 32.00"
kick = "0.05"
list_alpha = "0.50 1.00 3.00 10.00"
num_L = words(list_L)
num_V = words(list_V)
num_hz = words(list_hz)
num_alpha = words(list_alpha)

#set yrange [0.30:0.60]
set xrange [word(list_L,1)-0.3:word(list_L,num_L)+0.3]
set key out right box 3
unset colorbox
set size ratio 1.5 #0.8
#set ytics 0.3, 0.05, 0.6 format "%.2f"
set xtics word(list_L,1), 2, word(list_L,num_L) 
set xlabel "$L$"
set ylabel "$\\rangle log(Delta) \\langle"

i = 0 #alpha = 0.05
j = 0 # num_hz-1 #hz = 32.00
k = num_V-1 #V=3.00

#set table $DataSelected
#plot for [j=0:num_hz-1] file every :::i+num_alpha*j+num_alpha*num_hz*k::i+num_alpha*j+num_alpha*num_hz*k u 1:3:4:6:12:13 w table 
#unset table
#print "Table"
#print $DataSelected

set terminal pngcairo font ",13"
i = 0 #alpha = 0.05
j = num_hz-1 #hz = 32.00
k = num_V-1 #V=3.00
set output "figures/Scaling_Log_Gap_Difference_hz_L12.png"
set title "N_d = 2^{1-L/2}20480, T0 = ".T0.", J = ".J.", V = ".word(list_V,k+1).",\n kick = ".kick.", alpha = ".word(list_alpha,i+1)
plot for [j=0:num_hz-1] file every :::i+num_alpha*j+num_alpha*num_hz*k::i+num_alpha*j+num_alpha*num_hz*k u 1:12:13 w errorlines title "$hz = ".word(list_hz,j+1)."$" lc palette frac (j+0.0)/num_hz pt 4

set output "figures/Scaling_Log_Gap_Comparison_hz_L12.png"
plot for [j=0:num_hz-1] file every :::i+num_alpha*j+num_alpha*num_hz*k::i+num_alpha*j+num_alpha*num_hz*k u 1:14 w errorlines title "$hz = ".word(list_hz,j+1)."$, Delta_{pi}" lc palette frac (j+0.0)/num_hz pt 2 ,\
    for [j=0:num_hz-1] "" every :::i+num_alpha*j+num_alpha*num_hz*k::i+num_alpha*j+num_alpha*num_hz*k u 1:15 w errorlines title "$hz = ".word(list_hz,j+1)."$, Delta_0" lc palette frac (j+0.0)/num_hz pt 4 ,\
    log(pi/2**(x)) lw 1.2 lt 8 title "log(2*pi/2^L)"


i = 0 #alpha = 0.05
j = num_hz-1 #hz = 32.00
k = num_V-1 #V=3.00
set output "figures/Scaling_Log_Gap_Difference_alpha_L12.png"
set title "N_d = 2^{1-L/2}20480, T0 = ".T0.", J = ".J.", V = ".word(list_V,k+1).",\n kick = ".kick.", hz = ".word(list_hz,j+1)
plot for [i=0:num_alpha-1] file every :::i+num_alpha*j+num_alpha*num_hz*k::i+num_alpha*j+num_alpha*num_hz*k u 1:12:13 w errorlines title "$a = ".word(list_alpha,i+1)."$" lc palette frac (i+0.0)/num_alpha pt 4

set output "figures/Scaling_Log_Gap_Comparison_alpha_L12.png"
plot for [i=0:num_alpha-1] file every :::i+num_alpha*j+num_alpha*num_hz*k::i+num_alpha*j+num_alpha*num_hz*k u 1:14 w errorlines title "$a = ".word(list_alpha,i+1)."$, Delta_{pi}" lc palette frac (i+0.0)/num_alpha pt 2 ,\
     for [i=0:num_alpha-1] "" every :::i+num_alpha*j+num_alpha*num_hz*k::i+num_alpha*j+num_alpha*num_hz*k u 1:15 w errorlines title "$a = ".word(list_alpha,i+1)."$, Delta_0" lc palette frac (i+0.0)/num_alpha pt 4 ,\
    log(pi/2**(x)) lw 1.2 lt 8 title "log(2*pi/2^L)"


reset session

#################################################################################

set key out right box 3
unset colorbox

#y0 = 0.386
#y1 = 0.5295
#x0 = word(list_hz,1)
#x1 = word(list_hz,num_hz)
#set arrow from x0,y0 to x1,y0 nohead
#set arrow from x0,y1 to x1,y1 nohead
