##!/bin/bash
filestring="dynamics_swap_LR"

ulimit -s unlimited
make $filestring

n_threads=32
stack=6G
export OMP_NUM_THREADS=$n_threads 
export OMP_STACKSIZE=$stack

output="Swap_LR_dynamics_long.txt"
file_out="out.txt"
mv $file_out output/

###Choice of parameters
iterations_2=10240
steps="100000000"
list_nspin=$(seq 4 2 12)
list_J="0.1 0.01"
list_V="3.00"
list_hz="16.00"
list_kick="0.00"
list_alpha="0.50 3.00" #1.00 10.00
list_T0="1.00"

for nspin in $list_nspin
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  #iterations=10
  for kick in $list_kick
  do
    for T0 in $list_T0
    do
      for J in $list_J
      do
        for V in $list_V
        do
          for hz in $list_hz
          do
            for alpha in $list_alpha
            do
              cat > input.txt << *
$nspin
$iterations
$steps
$T0
$J
$V
$hz
$kick
$alpha
*
              echo ""
              echo "Input: "
              echo "nspin = $nspin  J = $J  V = $V  hz = $hz  epsilon = $kick  alpha = $alpha  T0 = $T0"
              echo "iterations = $iterations  steps = $steps  n_threads = $n_threads"
              echo ""
              echo "Executing file: "
  	          ./$filestring < input.txt | tee $file_out
            done
            echo ""
          done
        done
      done
    done
  done
done
