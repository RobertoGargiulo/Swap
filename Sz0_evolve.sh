##!/bin/bash
filestring="dense_Sz0"

make $filestring

#read nspin
#steps=200
iterations=80
n_threads=10
total_time=2000

export OMP_NUM_THREADS=$n_threads 
#j=0
for time_step in 0.50
do
  steps=`echo $total_time $time_step | awk '{print $1/$2}'`
  for nspin in 6 8 10 12 14 #16 #{2..18..2}
  do
    for kdim in 30 #50 100 #200
    do
      for J in 0.50
      do
        for V in 0.25 #1.5 2 2.5 3
        do
          for hz in $( seq 1.00 0.05 2.00 )
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
            echo "nspin = " $nspin
            echo "J = " $J " V = " $V " hz = " $hz
            echo "steps = " $steps "  time_step = " $time_step "  total_time = " $total_time
            echo "iterations = " $iterations "  kdim = " $kdim "  n_threads = " $n_threads
            avg=`cat data/magnetizations/Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin$nspin\_steps$steps\_time_step$time_step\_iterations$iterations\_J$J\_V$V\_hz$hz.txt | tail -3 | head -1`
            gap_ratio=`cat data/magnetizations/Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin$nspin\_steps$steps\_time_step$time_step\_iterations$iterations\_J$J\_V$V\_hz$hz.txt | tail -1`
            echo $J $V $hz $avg $gap_ratio >> data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin$nspin\_time_step$time_step\_steps$steps\_close_to_hz1.50.txt
          done
        done
      done
    done
  done
done
