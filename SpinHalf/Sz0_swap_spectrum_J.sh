##!/bin/bash
filestring="swap_spectrum_Sz0"

make $filestring

n_threads=30
export OMP_NUM_THREADS=$n_threads 
output="Swap_spectrum_J.txt"
echo "" >> $output
#mv $output ../Trash

iterations_2=5120
for nspin in 10 #4 6 8 10 #{6..16..2}
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  #iterations=1
  for kick in 0.00 #$(seq 0.00 0.05 0.20)
  do
    for period in 1.00
    do
      for J in 0.70 0.80 0.90 1.00 #$(seq 0.01 0.01 0.10; seq 0.20 0.10 1.00)
      do
        for V in 0.50 #$(seq 0.00 0.20 1.60)
        do
          for hz in 6.00 #0.01 $( seq 0.20 0.20 2.00; seq 4.00 4.00 16.00)
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
            file_out="out_J.txt"
  	        ./$filestring < input.txt | tee $file_out
            echo "nspin = $nspin "
            echo "J = $J  V = $V  hz = $hz  epsilon = $kick  period = $period"
            echo "iterations = $iterations   n_threads = $n_threads "
            gap_ratio=`grep -a -A1 "Ratio" $file_out | tail -1`
            log_gap=`grep -a -A1 "Logarithm" $file_out | tail -1`
            echo $nspin $J $V $hz $kick $period $gap_ratio $log_gap $iterations >> $output
          done
          echo "" >> $output
        done
      done
    done
  done
done
