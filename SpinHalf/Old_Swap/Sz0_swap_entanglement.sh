##!/bin/bash
filestring="swap_entanglement_Sz0"

make $filestring

#iterations=100 #1280
n_threads=16
export OMP_NUM_THREADS=$n_threads 
#j=0
#output="Swap_entanglement.txt"
#echo "" >> $output
#mv $output ../Trash

iterations_2=5120
for nspin in {8..14..2}
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  #iterations=100
  #iterations=320
  for kick in 0.00 #0.50
  do
    for period in 1.00
    do
      for J in 0.00 #1.00 #$(seq 0.00 0.05 0.50)
      do
        for V in 0.50 #$(seq 0.00 0.20 1.60)
        do
          for hz in 10.00 #0.00 2.00 $(seq 4.00 4.00 16.00)
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
            #entanglement=`grep -a -A1 "Information" $file_out | tail -1`
            #echo $nspin $J $V $hz $kick $period $entanglement $iterations >> $output
          done
          #echo "" >> $output
          echo ""
        done
      done
    done
  done
done
