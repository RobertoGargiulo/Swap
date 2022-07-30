##!/bin/bash
filestring="prova5"

make $filestring

#read nspin
steps=200
iterations=4
n_threads=2
total_time=100

export OMP_NUM_THREADS=$n_threads 
#j=0
for time_step in 0.5 0.05
do
  steps=`echo $total_time $time_step | awk '{print $1/$2}'`
  echo "steps = " $steps
  for nspin in 12 #{2..18..2}
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
