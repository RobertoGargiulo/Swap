##!/bin/bash
filestring="prova3"

make $filestring

#export OMP_NUM_THREADS=8
#read ncores
#read nspin
#ncores=4
nspin=12
steps=10000
iterations=100
time_step=0.5
#for nspin in 2 4 6 8
#do
kdim=100
for n_threads in 1 2 4 8
do
  export OMP_NUM_THREADS=$n_threads 
  for kick in 0.00 #0.05 0.10 0.20 0.40
  do
    for J in 0.5 #1 2 #1 1.5 2
    do
      for V in 0.5 #1 #0 0.5 1 1.5 2 3
      do
        for hz in 3 #0 1 2.5 3 4 6
        do
          cat > input.txt << *
$nspin
$iterations
$steps
$kdim
$time_step
$J
$V
$hz
$kick
*
	        ./$filestring < input.txt
          #./mag_avg < input.txt
        done
      done
    done
  done
done
