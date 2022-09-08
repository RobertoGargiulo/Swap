##!/bin/bash
filestring="swap_dense_Sz0"

make $filestring

#iterations=100 #1280
n_threads=10
total_time=10000
export OMP_NUM_THREADS=$n_threads 
#j=0
#output="Swap_varying_V.txt"
#mv $output ../Trash

iterations_2=2560
for nspin in {2..12..2}
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  for kick in 0.00
  do
    for period in 1.00
    do
      steps=`echo $total_time $period | awk '{print $1/$2}'`
      for J in 0.05 #0.50 1.57 3.14
      do
        for V in 0.50 #$(seq 0.00 0.30 3.00; seq 3.05 0.05 3.15)
        do
          for hz in 6.00 #$(seq 0.00 1.00 6.00)
          do
            cat > input.txt << *
$nspin
$iterations
$steps
$period
$J
$V
$hz
$kick
*
            file_out="out.txt"
  	        ./$filestring < input.txt | tee $file_out
            echo "nspin = $nspin "
            echo "J = $J  V = $V  hz = $hz  epsilon = $kick"
            echo "steps = $steps   period = $period   total_time = $total_time "
            echo "iterations = $iterations   n_threads = $n_threads "
            avg=`grep -a -A1 "Averages" $file_out | tail -1`
            #echo $nspin $J $V $hz $kick $period $avg $gap_ratio $iterations $steps >> $output
          done
          #echo "" >> $output
        done
      done
    done
  done
done
