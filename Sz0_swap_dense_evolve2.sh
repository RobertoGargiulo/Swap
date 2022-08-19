##!/bin/bash
filestring="swap_dense_Sz0"

make $filestring

#iterations=100 #1280
n_threads=20
total_time=1000000
export OMP_NUM_THREADS=$n_threads 
output="Swap_decay_times_varying_eps.txt"

iterations_2=5120
for nspin in 2 4 6 8
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  for kick in $(seq 0.00 0.05 0.20) 
  do
    for period in 1.00 
    do
      steps=`echo $total_time $period | awk '{print $1/$2}'`
      for J in 0.05 #$(seq 0.00 0.04 0.24) 
      do
        for V in 0.50 #$(seq 0.00 0.40 1.60) 
        do
          for hz in 12.00 #0.01 2.00 $(seq 4.00 4.00 16.00)
          do
            cat > input.txt << *
$nspin
$iterations
$steps
$period
$J
$V
$hz
$kick
*
            file_out="out5.txt"
  	        ./$filestring < input.txt | tee $file_out
            echo "nspin = $nspin "
            echo "J = $J  V = $V  hz = $hz  epsilon = $kick"
            echo "period = $period   steps = $steps   total_time = $total_time "
            echo "iterations = $iterations   n_threads = $n_threads "
            avg=`grep -a -A1 "Imbalance" $file_out | tail -1`
            #gap_ratio=`grep -a -A1 "Gap Ratio" $file_out | tail -1`
            decay_times=`grep -a -A1 "Decay Time" $file_out | tail -1`
            echo $nspin $J $V $hz $kick $period $avg $decay_times $iterations $steps >> $output
          done
          echo "" >> $output
        done
      done
    done
  done
done
