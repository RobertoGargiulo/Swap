##!/bin/bash
filestring="swap_dense_Sz0"

make $filestring

#read nspin
#steps=200
iterations=20 #1280
n_threads=2
#total_time=2000
#kick=0
steps=8000
export OMP_NUM_THREADS=$n_threads 
#j=0
  for nspin in 8 #10 12 #{6..16..2}
  do
for kick in 0 #0.1 0.5
do
for time_step in 1 #0.1 1 10.00
do
  #steps=`echo $total_time $time_step | awk '{print $1/$2}'`
    #iterations=`echo $iterations | awk '{print $1/2}'`
    for J in 0.01 #0.1  #0.50
    do
      for V in 0.25 #6 10 #3 5 #1.5 2 2.5 3
      do
        #output="data/phases/PT_Sz0_DENSE_SWAP_hz_Disorder_AVG_FLUCT_Imbalance_nspin${nspin}_time_step${time_step}_steps${steps}_iterations${iterations}_J${J}_V${V}_up_to_hz6.txt"
        output="out_Swap_trial3.txt"
        for hz in 4.5 #6 12.00 #$( seq 0.50 0.50 6.00) #; seq 4.10 0.10 6.00 )
        do
          cat > input.txt << *
$nspin
$iterations
$steps
$time_step
$J
$V
$hz
$kick
*
  	      ./$filestring < input.txt | tee out3.txt 
          input="data/magnetizations/Sz0_DENSE_SWAP_hz_Disorder_AVG_FLUCT_Imbalance_nspin${nspin}_steps${steps}_time_step${time_step}_iterations${iterations}_J${J}_V${V}_hz${hz}.txt"
          echo "nspin = $nspin "
          echo "J = $J  V = $V  hz = $hz  epsilon = $kick"
          echo "steps = $steps   period = $time_step   total_time = $total_time "
          echo "iterations = $iterations   n_threads = $n_threads "
          avg=`grep -A1 "Averages" out3.txt | tail -1`
          gap_ratio=`grep -A1 "Ratio" out3.txt | tail -1`
          echo $nspin $J $V $hz $kick $time_step $avg $gap_ratio >> $output
        done
      done
    done
  done
done
done
