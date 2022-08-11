##!/bin/bash
filestring="swap_dense_Sz0"

make $filestring

#iterations=100 #1280
n_threads=10
total_time=4000
export OMP_NUM_THREADS=$n_threads 
#j=0
output="out_Swap_PD_data.txt"
#mv $output ../Trash

iterations_2=2560
for nspin in 12 #2 4 6 8 10 #{6..16..2}
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  for kick in $(seq 0.00 0.05 0.20)
  do
    for period in 1 #0.1 1 10.00
    do
      steps=`echo $total_time $period | awk '{print $1/$2}'`
      for J in $(seq 0.00 0.10 0.50)  #0.50
      do
        for V in $(seq 0.00 0.50 3.00) #6 10 #3 5 #1.5 2 2.5 3
        do
          #output="data/phases/PT_Sz0_DENSE_SWAP_hz_Disorder_AVG_FLUCT_Imbalance_nspin${nspin}_period${period}_steps${steps}_iterations${iterations}_J${J}_V${V}_up_to_hz6.txt"
          for hz in $(seq 0.00 1.00 6.00) #$( seq 0.50 0.50 6.00) #; seq 4.10 0.10 6.00 )
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
          file_out="out.txt"
  	      ./$filestring < input.txt | tee $file_out
          echo "nspin = $nspin "
          echo "J = $J  V = $V  hz = $hz  epsilon = $kick"
          echo "steps = $steps   period = $period   total_time = $total_time "
          echo "iterations = $iterations   n_threads = $n_threads "
          avg=`grep -A1 "Averages" $file_out | tail -1`
          gap_ratio=`grep -A1 "Ratio" $file_out | tail -1`
          echo $nspin $J $V $hz $kick $period $avg $gap_ratio $iterations $steps >> $output
        done
      done
    done
  done
done
done
