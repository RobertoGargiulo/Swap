##!/bin/bash
filestring="dense_Sz0"

make $filestring

#read nspin
#steps=200
iterations=1280
n_threads=10
total_time=2000

export OMP_NUM_THREADS=$n_threads 
#j=0
for time_step in 0.50
do
  steps=`echo $total_time $time_step | awk '{print $1/$2}'`
  for J in 0.50
  do
    for V in 0.25 #1.5 2 2.5 3
    do
      for hz in 0.05 1.00 3.50 6.00 10.00 30.00  #$( seq 11.0 1.00 20.00 )
      do
        iterations=1280
        output="data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_time_step${time_step}_steps${steps}_J${J}_V${V}_hz${hz}_up_to_nspin12.txt"
        for nspin in 14 #{6..12..2}
        do
          #iterations=`echo $iterations | awk '{print $1/2}'`
          iterations=40
          cat > input.txt << *
$nspin
$iterations
$steps
$time_step
$J
$V
$hz
*
  	      ./$filestring < input.txt | tee out2.txt 
          input="data/magnetizations/Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin${nspin}_steps${steps}_time_step${time_step}_iterations${iterations}_J${J}_V${V}_hz${hz}.txt"
          echo "nspin = $nspin "
          echo "J = $J  V = $V  hz = $hz "
          echo "steps = $steps   time_step = $time_step   total_time = $total_time "
          echo "iterations = $iterations   n_threads = $n_threads "
          avg=`grep -A1 "Averages" out2.txt | tail -1`
          gap_ratio=`grep -A1 "Ratio" out2.txt | tail -1`
          echo $nspin $J $V $avg $gap_ratio $iterations >> $output
        done
      done
    done
  done
done
