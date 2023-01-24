
set terminal png # tikz color # standalone
set size ratio 0.5

file = "data.txt"
set output "figures/figure.png"
set yrange [-0.1:1.1]
set xrange [-0.1:0.5]
plot file u 1:3 ps 4
