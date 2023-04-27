##!/bin/bash
filestring="test_swap_LRZZ"

make $filestring

n_threads=1
export OMP_NUM_THREADS=$n_threads 
#output="Swap_spectrum_ergodic_kick_density0.05.txt"
#output="Swap_spectrum_check_GOE_for_entanglement.txt"
output="Swap_LR_degeneracies_solvable_point_check.txt"
echo "" >> $output
#mv $output ../Trash

iterations_2=5120 #20480
for nspin in 12 #{6..14..2}
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  iterations=10
  for kick in 0.00 #-0.85 -0.79 $(seq -0.75 0.05 0.75) 0.79 0.85
  do
    for period in 1.00
    do
      for J in 0.00 #1.00 #$(seq 0.00 0.05 0.30; seq 0.40 0.20 1.00)
      do
        for V in 1.00 #$(seq 0.00 0.20 1.60)
        do
          #for idx in {0..10..1}
          for hz in 4.00
          do
            #hz=`echo $(printf "%.2f" $(echo "0.04 * 2^($idx)" | bc) )`
            for alpha in 0.50 1.00 3.00 5.00 10.00
            do
              cat > input.txt << *
$nspin
$iterations
$period
$J
$V
$hz
$kick
$alpha
*
              file_out="out.txt"
  	          ./$filestring < input.txt | tee $file_out
              echo "nspin = $nspin "
              echo "J = $J  V = $V  hz = $hz  epsilon = $kick  period = $period  alpha = $alpha"
              echo "iterations = $iterations   n_threads = $n_threads "
              gap_ratio=`grep -a -A2 "Ratio" $file_out | tail -1`
              log_gap=`grep -a -A2 "Logarithm" $file_out | tail -1`
              echo $nspin $J $V $hz $kick $alpha $period $gap_ratio $log_gap $iterations >> $output
            done
            echo ""
          done
          #echo "" >> $output
        done
      done
    done
  done
  echo "" >> $output
done
