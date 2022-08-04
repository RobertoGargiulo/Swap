##!/bin/bash
filestring="spectrum_Sz0"

make $filestring

#read nspin
#steps=200
iterations=80
n_threads=10

export OMP_NUM_THREADS=$n_threads 
#j=0
  for nspin in 6 8 10 12 #6 8 10 12 14 #16 #{2..18..2}
  do
    for J in 0.50
    do
      for V in $( seq 0.00 0.05 3.00 ) #1.5 2 2.5 3
      do
        for hz in $( seq 0.00 0.05 6.00 )
        do
          cat > input.txt << *
$nspin
$iterations
$J
$V
$hz
*
  	      ./$filestring < input.txt | tee out.txt
          echo "nspin = " $nspin
          echo "J = " $J " V = " $V " hz = " $hz
          echo "iterations = " $iterations "  n_threads = " $n_threads
          gap_ratio=`cat out.txt | grep -A1 Ratio | tail -1`
          echo $J $V $hz $gap_ratio >> data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_Gap_Ratio_nspin$nspin\_iterations$iterations\_J$J\_up_to_V3_hz6.txt
          done
        done
      done
    done
