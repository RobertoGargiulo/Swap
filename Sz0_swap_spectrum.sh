##!/bin/bash
filestring="swap_spectrum_Sz0"

make $filestring

#iterations=100 #1280
n_threads=10
export OMP_NUM_THREADS=$n_threads 
#j=0
output="Swap_entanglement_trial.txt"
#mv $output ../Trash
pi=3.14
pi2=1.57
endV=`echo $pi | awk '{print 4*$1}'`

iterations_2=320
for nspin in 2 4 6 8 #2 4 6 8 10 #{6..16..2}
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  for kick in $(seq 0.00 0.05 0.20)
  do
    for period in 1.00
    do
      for J in 0.05 0.50 #1.57 3.14
      do
        for V in 0.25 0.50 #$(seq 0.00 0.40 1.60)
        do
          #output="data/phases/PT_Sz0_DENSE_SWAP_hz_Disorder_AVG_FLUCT_Imbalance_nspin${nspin}_period${period}_steps${steps}_iterations${iterations}_J${J}_V${V}_up_to_hz6.txt"
          for hz in 3.00 12.00 #$(seq 0.00 4.00 16.00)
          do
            cat > input.txt << *
$nspin
$iterations
$period
$J
$V
$hz
$kick
*
            file_out="out.txt"
  	        ./$filestring < input.txt | tee $file_out
            echo "nspin = $nspin "
            echo "J = $J  V = $V  hz = $hz  epsilon = $kick  period = $period"
            echo "iterations = $iterations   n_threads = $n_threads "
            gap_ratio=`grep -a -A1 "Ratio" $file_out | tail -1`
            entanglement=`grep -a -A1 "Information" $file_out | tail -1`
            echo $nspin $J $V $hz $kick $period $gap_ratio $entanglement $iterations >> $output
          done
          echo "" >> $output
        done
      done
    done
  done
done
