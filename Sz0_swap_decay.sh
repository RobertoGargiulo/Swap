##!/bin/bash
filestring="swap_decay"

make $filestring

n_threads=60
#total_time=1000000
export OMP_NUM_THREADS=$n_threads 
output="Swap_decay_times_Large_IMB_LI.txt"
#mv $output ../Trash

iterations_2=5120
for nspin in 4 6 8 10 12
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  #iterations=100
  for kick in $(seq 0.00 0.05 0.20) 
  do
    for period in 1.00 
    do
      for J in $(seq 0.05 0.05 0.30) 
      do
        for V in $(seq 0.00 0.40 1.60) 
        do
          for hz in 0.01 2.00 $(seq 5.00 5.00 20.00)
          do
            total_time=$((5*10**(nspin/2)))
            n_periods=10
            steps=`echo $total_time $period $n_periods | awk '{print $1/$2/$3}'`

            cat > input.txt << *
$nspin
$iterations
$steps
$n_periods
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
            echo "period = $period  n_periods = $n_periods  steps = $steps   total_time = $total_time "
            echo "iterations = $iterations   n_threads = $n_threads "
            decay_times=`grep -a -A1 "Decay Time" $file_out | tail -1`
            #echo $nspin $J $V $hz $kick $period $decay_times $iterations $total_time $n_periods $steps >> $output
          done
          #echo "" >> $output
          echo ""
        done
      done
    done
  done
done
