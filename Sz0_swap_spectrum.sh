##!/bin/bash
filestring="swap_spectrum_Sz0"

make $filestring

#read nspin
#steps=200
iterations=1280
n_threads=10

export OMP_NUM_THREADS=$n_threads 
#j=0
for nspin in 6 8 10 12 # 14 16 #{2..18..2}
do
  iterations=`echo $iterations | awk '{print $1/2}'`
  for kick in 0 0.1 0.5
  do
    for time_step in 0.10 1.00 10.00
    do
      for J in 0.50
      do
        for V in $( seq 0.00 0.20 6.00 )
        do
          for hz in $( seq 0.00 0.40 12.00 )
          do
            output="data/phases/PT_Sz0_DENSE_SWAP_hz_Disorder_AVG_Gap_Ratio_nspin${nspin}_iterations${iterations}_period${time_step}_kick${kick}_J${J}_up_to_V5_hz10.txt"
            cat > input.txt << *
$nspin
$iterations
$time_step
$J
$V
$hz
$kick
*
  	        ./$filestring < input.txt | tee out.txt
            echo "nspin = $nspin"
            echo "J = $J  V = $V  hz = $hz  epsilon = $kick"
            echo "period = $time_step"
            echo "iterations = $iterations   n_threads = $n_threads "
            gap_ratio=`cat out.txt | grep -A1 Ratio | tail -1`
            echo $J $V $hz  $gap_ratio >> $output
            done
          done
        done
      done
    done
    echo "" >> $output
  done
