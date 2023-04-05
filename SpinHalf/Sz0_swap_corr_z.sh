##!/bin/bash
filestring="swap_corr_z_Sz0"

make $filestring

#iterations=100 #1280
n_threads=8
export OMP_NUM_THREADS=$n_threads 
#j=0
output="Swap_correlations.txt"
echo "----------------------" >> $output
#echo "" >> $output
#mv $output ../Trash

iterations_2=5120
for nspin in {2..12..2}
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  #iterations=1
  for kick in 0.05 #0.50
  do
    for period in 1.00
    do
      for J in 0.10 #$(seq 0.00 0.10 0.50) #$(seq 0.00 0.05 0.50)
      do
        for V in 0.50 #$(seq 0.00 0.40 1.60) #$(seq 0.00 0.20 1.60)
        do
          for hz in 0.00 0.01 0.10 0.20 0.50 1.00 2.00 $(seq 4.00 4.00 20.00) #$(seq 4.00 4.00 16.00)
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
            correlation=`grep -a -A1 "average correlation" $file_out | tail -1`
            norm_correlation=`grep -a -A1 "normalized" $file_out | tail -1`
            echo $nspin $J $V $hz $kick $period $correlation $norm_correlation $iterations >> $output
          done
          echo "" >> $output
          echo ""
        done
      done
    done
  done
done
