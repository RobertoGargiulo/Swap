##!/bin/bash
filestring="prova3"

make $filestring

#export OMP_NUM_THREADS=8
#read ncores
#read nspin
#ncores=4
#nspin=8
steps=4000
iterations=30
#time_step=0.5
n_threads=3
#j=0
export OMP_NUM_THREADS=$n_threads 
#kdim=50
for time_step in 0.5 #0.05 0.5 1 2
do
for nspin in 10 #{2..18..2}
do
  for kdim in 30 #50 100 200
  do
    for J in 0.5 #1 2
    do
      for V in 1 #1 #0 0.5 1 1.5 2 2.5 3
      do
        for hz in 3 #0 1 2 3 4 5 6
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
          #j=`echo $j | awk '{print $1 + 1}'`
          #echo "loop " $j
	        ./$filestring < input.txt
          #cat input.txt
        done
      done
    done
  done
done
done
