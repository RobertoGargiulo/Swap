##!/bin/bash
filestring="prova2"

make $filestring

#export OMP_NUM_THREADS=8
#read ncores
#read nspin
#ncores=4
nspin=10
steps=4000
iterations=30
time_step=0.5
#for nspin in 2 4 6 8
#do
kdim=100
n_threads=2
export OMP_NUM_THREADS=$n_threads 
    for J in 0.5 #1 2
    do
      for V in 1 #0 0.13 0.38 0.5 1 2
      do
        for hz in 3 #0 1.25 1.5 2 2.5 3 4 6
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
*
	        ./$filestring < input.txt
        done
      done
    done
