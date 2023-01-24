##!/bin/bash
filestring="swap_overlap"

make $filestring

#iterations=100 #1280
n_threads=10
export OMP_NUM_THREADS=$n_threads 
#j=0
#output="Swap_overlap_Neel.txt"
#echo "" >> $output
#mv $output ../Trash

iterations_2=20480
for nspin in {2..14..2}
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  #iterations=100
  for kick in 0.05 #$(seq 0.00 0.05 0.20)
  do
    for period in 1.00
    do
      for J in $(seq 0.05 0.05 0.30) 0.40 0.50 1.00 2.00 
      do
        for V in 0.50 #$(seq 0.00 0.20 1.60)
        do
          for hz in 6.00 #0.01 2.00 $(seq 4.00 4.00 16.00)
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
            file_out="out3.txt"
  	        ./$filestring < input.txt | tee $file_out
            echo "nspin = $nspin "
            echo "J = $J  V = $V  hz = $hz  epsilon = $kick  period = $period"
            echo "iterations = $iterations   n_threads = $n_threads "
            #entanglement=`grep -a -A1 "Information" $file_out | tail -1`
            #echo $nspin $J $V $hz $kick $period $entanglement $iterations >> $output
            echo ""
          done
          #echo "" >> $output
        done
      done
    done
  done
done
