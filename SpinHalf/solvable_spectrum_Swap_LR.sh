##!/bin/bash
filestring="solvable_spectrum_swap_LR"

make $filestring

n_threads=1
export OMP_NUM_THREADS=$n_threads 

output="Solvable_Swap_LR_spectrum.txt"
mv $output data/eigen/
file_out="out2.txt"
file_sort="out_sort2.txt"
mv temp $file_out $file_sort raw_$output output/

block=0
iterations_2=5120 #20480
for nspin in {4..18..2}
do
  block=$(( $block + 1 ))
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  #iterations=1
  for period in 1.00
  do
    for V in 2.00 #$(seq 0.00 0.20 1.60)
    do
      #for idx in {0..10..1}
      for hz in 1.00 2.00 4.00 8.00 16.00
      do
        #hz=`echo $(printf "%.2f" $(echo "0.04 * 2^($idx)" | bc) )`
        for alpha in 0.50 1.00 3.00 $(seq 4.00 1.00 10.00)
        do
          cat > input.txt << *
$nspin
$iterations
$period
$V
$hz
$alpha
*
          ./$filestring < input.txt | tee $file_out
          echo "nspin = $nspin"
          echo "V = $V  hz = $hz  period = $period  alpha = $alpha"
          echo "iterations = $iterations   n_threads = $n_threads"
          gap_ratio=`grep -a -A2 "Ratio" $file_out | tail -1`
          log_gap=`grep -a -A2 "Logarithm" $file_out | tail -1`
          ord_gap=`grep -a -A2 "Ordinary" $file_out | tail -1`
          tot_deg=`grep -a -A1 "degenerate" $file_out | tail -1`
          echo $nspin $V $hz $alpha $period $gap_ratio $log_gap $ord_gap $tot_deg $iterations >> raw_$output

        done
        echo "" >> raw_$output
        echo ""
      done
    done
  done
  column -t < raw_$output >> temp
  nparam=3
  cat > input_sort.txt << *
$nparam
1
3
4
4
*
  ./sort_data.sh temp < input_sort.txt | tee $file_sort
  mv sort_col4_temp $output
  mv temp output/

  cat > input_sort.txt << *
$nparam
3
4
1
$block
*
  ./sort_data.sh $output < input_sort.txt | tee $file_sort
done

mv temp $file_out $file_sort raw_$output output/
