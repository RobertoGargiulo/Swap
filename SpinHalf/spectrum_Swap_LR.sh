##!/bin/bash
filestring="spectrum_swap_LR"

make $filestring

n_threads=10
export OMP_NUM_THREADS=$n_threads 

output="Swap_LR_spectrum_L12_V_hz_alpha_large_N_dis.txt"
mv $output data/eigen/
file_out="out.txt"
file_sort="out_sort.txt"
mv temp $file_out $file_sort raw_$output output/

block=0
iterations_2=20480
for nspin in {4..12..2}
do
  block=$(( $block + 1 ))
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  #iterations=10
  for kick in 0.05 #-0.85 -0.79 $(seq -0.75 0.05 0.75) 0.79 0.85
  do
    for period in 1.00
    do
      for J in 0.10 #1.00 #$(seq 0.00 0.05 0.30; seq 0.40 0.20 1.00)
      do
        for V in 1.00 2.00 3.00 #$(seq 0.00 0.20 1.60)
        do
          #for idx in {0..10..1}
          for hz in 2.00 4.00 8.00 16.00 32.00
          do
            #hz=`echo $(printf "%.2f" $(echo "0.04 * 2^($idx)" | bc) )`
            for alpha in 0.50 1.00 3.00 10.00 #20.00
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
              echo "nspin = $nspin"
              echo "J = $J  V = $V  hz = $hz  epsilon = $kick  period = $period  alpha = $alpha"
              echo "iterations = $iterations   n_threads = $n_threads"
              gap_ratio=`grep -a -A2 "Ratio" $file_out | tail -1`
              log_gap=`grep -a -A2 "Logarithm" $file_out | tail -1`
              echo $nspin $J $V $hz $kick $alpha $period $gap_ratio $log_gap $iterations >> raw_$output
            done
            echo ""
            echo "" >> raw_$output
          done
          #echo "" >> $output
        done
      done
    done
  done
  column -t < raw_$output >> temp
  nparam=4
  cat > input_sort.txt << *
$nparam
1
3
4
6
4
*
  ./sort_data.sh temp < input_sort.txt | tee $file_sort
  mv sort_col6_temp $output
  mv temp output/

  cat > input_sort.txt << *
$nparam
3
4
6
1
$block
*
  ./sort_data.sh $output < input_sort.txt | tee $file_sort
done

mv temp $file_out $file_sort raw_$output output/
