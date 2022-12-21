##!/bin/bash
filestring="swap_spectrum_Sz0"

make $filestring

n_threads=20
export OMP_NUM_THREADS=$n_threads 
#output="Swap_spectrum_ergodic_kick_density0.05.txt"
#output="Swap_spectrum_check_GOE_for_entanglement.txt"
output="Swap_solvable_point_check.txt"
echo "" >> $output
#mv $output ../Trash

iterations_2=5120 #20480
for nspin in {8..14..2}
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  #iterations=1
  for kick in 0.00 #-0.85 -0.79 $(seq -0.75 0.05 0.75) 0.79 0.85
  do
    for period in 1.00
    do
      for J in 0.00 #1.00 #$(seq 0.00 0.05 0.30; seq 0.40 0.20 1.00)
      do
        for V in 0.50 #1.00 #$(seq 0.00 0.20 1.60)
        do
          #for idx in {0..10..1}
          for hz in 10.00 #4.00
          do
            #hz=`echo $(printf "%.2f" $(echo "0.04 * 2^($idx)" | bc) )`
            
            cat > input.txt << *
$nspin
$iterations
$period
$J
$V
$hz
$kick
*
            file_out="out2.txt"
  	        ./$filestring < input.txt | tee $file_out
            echo "nspin = $nspin "
            echo "J = $J  V = $V  hz = $hz  epsilon = $kick  period = $period"
            echo "iterations = $iterations   n_threads = $n_threads "
            gap_ratio=`grep -a -A1 "Ratio" $file_out | tail -1`
            log_gap=`grep -a -A1 "Logarithm" $file_out | tail -1`
            echo $nspin $J $V $hz $kick $period $gap_ratio $log_gap $iterations >> $output
          done
          echo ""
          #echo "" >> $output
        done
      done
    done
  done
  echo "" >> $output
done
