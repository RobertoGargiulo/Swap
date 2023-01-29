##!/bin/bash
filestring="dynamics_MBL"

make $filestring

n_threads=10
#total_time=1000000
export OMP_NUM_THREADS=$n_threads 


iterations_2=1280 #5120
#echo "nspin = " ; read nspin
period=1.00
for nspin in {2..8..2}
do
  #iterations=$iterations_2
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  for kick in 0.00
  do
    for J in $(seq 0.20 0.20 1.00; seq 1.20 0.40 2.00)
    do
      for V in $(seq 0.00 0.30 3.00) 
      do
        for hz in $(seq 0.00 2.00 10.00)
        do
          for steps in 10000
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
            file_out="out_MBL_grid_J_V_hz.txt"
  	        ./$filestring < input.txt | tee $file_out
            echo "nspin = $nspin"
            echo "J = $J  V = $V  hz = $hz  epsilon = $kick"
            echo "period = $period  steps = $steps" 
            echo "iterations = $iterations  n_threads = $n_threads "
            #echo $nspin $J $V $hz $kick $period $decay_times $iterations $total_time $n_periods $steps >> $output
          done
          #echo "" >> $output
          echo ""
        done
      done
    done
  done
done
