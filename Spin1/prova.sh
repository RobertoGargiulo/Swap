##!/bin/bash
filestring="prova"

make $filestring

n_threads=1
#total_time=1000000
export OMP_NUM_THREADS=$n_threads 


iterations_2=5120
echo "nspin = " ; read nspin
  iterations=1
  for kick in 0.05
  do
    for period in 1.00 
    do
      for J in 0.05
      do
        for V in 0.50
        do
          for hz in 10.00
          do
            for steps in 1
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
  	        ./$filestring < input.txt
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
#done
