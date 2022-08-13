##!/bin/bash
filestring="swap_dense_Sz0"

make $filestring

#iterations=100 #1280
n_threads=30
total_time=100000
export OMP_NUM_THREADS=$n_threads 
output="Swap_decay_times.txt"

iterations_2=5120
for nspin in 2 4 6 8 10
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  for kick in 0.00 
  do
    for period in 1.00 
    do
      steps=`echo $total_time $period | awk '{print $1/$2}'`
      for J in $(seq 0.02 0.02 0.10) 
      do
        for V in 0.25 0.50 #$(seq 0.00 0.10 4.00) 
        do
          for hz in 3.00 6.00 12.00
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
            file_out="out2.txt"
  	        ./$filestring < input.txt | tee $file_out
            echo "nspin = $nspin "
            echo "J = $J  V = $V  hz = $hz  epsilon = $kick"
            echo "steps = $steps   period = $period   total_time = $total_time "
            echo "iterations = $iterations   n_threads = $n_threads "
            avg=`grep -a -A1 "Imbalance" $file_out | tail -1`
            gap_ratio=`grep -a -A1 "Gap Ratio" $file_out | tail -1`
            decay_times=`grep -a -A1 "Decay Time" $file_out | tail -1`
            echo $nspin $J $V $hz $kick $period $avg $decay_times $iterations $steps >> $output
          done
          echo "" >> $output
        done
        echo "" >> $output
      done
      echo "" >> $output
    done
    echo "" >> $output
  done
  echo "" >> $output
done
