########## Log Pi-Gap Swap Small J wrt J #############

file = "Swap_LR_decay_kick0_L12.txt"

list_L = "4 6 8 10 12"
T0 = "1.00"
V = "3.00"
hz = "16.00"
list_J = "0.0625 0.125 0.25 0.50 1.0"
kick = "0.0"
list_alpha = "0.50 3.00"
num_L = words(list_L)
num_J = words(list_J)
num_alpha = words(list_alpha)
f(x,y) = x+num_alpha*y ### Goes through the blocks with varying J, where x,y change alpha, J respectively

do for [j=0:num_alpha-1] {
print j
alpha = word(list_alpha,j+1) 
print "alpha = ", alpha
set table $DataSelected
plot for [i=0:num_J-1] file every :::f(j,i)::f(j,i) u 1:2:3:4:5:6:7:8:11:($7*$11) w table
unset table
print "Table"
print $DataSelected
}

#set xrange [word(list_L,1)-0.3:word(list_L,num_L)+0.3]
#set xrange [word(list_lambda,1):word(list_L,num_L)+0.3]
set key out right box 3
unset colorbox
set size ratio 0.8
set logscale y
set xtics 4, 2, 12
set xlabel "$L$"
set ylabel "$\\log_{10}(\\Delta)$"


do for [j=0:num_alpha-1] {

print j
alpha = word(list_alpha,j+1) 
print "alpha = ", alpha

set margin 0
set key width 1.1 spacing 1.1


fmt = ".tex"
#set terminal pngcairo dashed font ",13"
set terminal tikz color standalone scale 0.65, 0.65

set output "figures/Swap_LR_Decay_Times_wrt_J_kick0_alpha".alpha."_L".word(list_L,num_L)."".fmt
#set title "N_d = 2^{1-L/2}20480, T0 = ".T0.", V = ".V.", kick = ".kick.",\n hz = ".hz.", alpha = ".alpha
set notitle
plot for [i=0:num_J-1] file every :::f(j,i)::f(j,i) u 1:($7*$11):($8*$11) w errorlines title "$J = ".word(list_J,i+1)."$" lc palette frac (i+0.0)/num_L pt 4

}

reset session
