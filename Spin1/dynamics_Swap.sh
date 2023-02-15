##!/bin/bash
filestring="dynamics_Swap"

make $filestring

n_threads=20
#total_time=1000000
export OMP_NUM_THREADS=$n_threads 


iterations_2=640 #2560 #5120
#echo "nspin = " ; read nspin
period=1.00
for nspin in {2..8..2}
do
  #iterations=10
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  for kick in 0.05 #0.05 0.10 #0.20
  do
    for J in 0.05 #0.05 0.10 0.20 #0.00 0.25 0.50
    do
      for V in 3.00 #$(seq 0.50 0.50 5.00)
      do
        for hz in 5.00
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
            file_out="out_Swap.txt"
  	        ./$filestring < input.txt | tee $file_out
            echo "nspin = $nspin"
            echo "J = $J  V = $V  hz = $hz  epsilon = $kick"
            echo "period = $period  steps = $steps" 
            echo "iterations = $iterations  n_threads = $n_threads"
            #echo $nspin $J $V $hz $kick $period $decay_times $iterations $total_time $n_periods $steps >> $output
          done
          #echo "" >> $output
          echo ""
        done
      done
    done
  done
done

# list_kick = 0.00 0.05 0.10 0.20 (J0.05, V10.00, hz20.00)
# list_J = 0.00 0.05 0.10 0.20 (kick0.05, V10.00, hz20.00)
# list_V = 5.00 10.00 15.00 20.00 (kick0.05, J0.05, hz20.00)
# list_hz = 5.00 10.00 20.00 30.00 (kick0.05, J0.05, V10.00)
# list_L = 2 4 6 8 (mostly 8)
