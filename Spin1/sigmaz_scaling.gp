reset session

set size ratio 0.5
set margins -1,2,0,0

Data = "data/dynamics/sigmaz_MBL_Neel_nspin8_steps10000_period1.00_n_disorder160_J0.60_V2.70_hz4.00.txt"

set terminal png
set output "figures/figure.png"

n_dis2 = 1280

L = "8"
steps = "10000"
T0 = 1
n_dis = 2

set title Data

set logscale x
#set yrange [0:1]
plot Data u 1:( (sum[i=1:int(L)] (-1)**i * column(i+1))/L ) w lp lw 1.6 title "MBL" ,\



