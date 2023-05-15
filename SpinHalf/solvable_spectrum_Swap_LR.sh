##!/bin/bash
filestring="solvable_spectrum_swap_LR"

make $filestring

n_threads=4
export OMP_NUM_THREADS=$n_threads 

output="Solvable_Swap_LR_spectrum.txt"
mv $output data/eigen/
file_out="out2.txt"
file_sort="out_sort2.txt"
mv temp $file_out $file_sort raw_$output output/

###Choice of parameters
iterations_2=20480
list_nspin=$(seq 4 2 18)
list_V="3.00"
list_hz="16.00"
list_alpha="0.50 1.00 3.00 10.00 50.00" #20.00
list_T0="1.00"


nums=()
nums+=(`wc -w <<< $list_nspin`)
nums+=(`wc -w <<< $list_V`)
nums+=(`wc -w <<< $list_hz`)
nums+=(`wc -w <<< $list_alpha`)
nums+=(`wc -w <<< $list_T0`)

cols=`echo "${nums[@]}" | awk '{ indexes=""; for(i=1;i<=NF;i++) {if($i>1) {indexes=indexes i " "}} print indexes }'`
cols=($cols)
nparam=${#cols[@]}
echo "nums = ${nums[@]}"
echo "cols = "${cols[@]}
echo "nparam = "$nparam

block=0
for nspin in $list_nspin
do
  block=$(( $block + 1 ))
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  #iterations=1
  for T0 in $list_T0
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
$T0
$V
$hz
$alpha
*
          ./$filestring < input.txt | tee $file_out
          echo "nspin = $nspin"
          echo "V = $V  hz = $hz  T0 = $T0  alpha = $alpha"
          echo "iterations = $iterations   n_threads = $n_threads"
          gap_ratio=`grep -a -A2 "Ratio" $file_out | tail -1`
          log_gap=`grep -a -A2 "pi-Logarithmic" $file_out | tail -1`
          shift_log_gap=`grep -a -A2 "Half Shifted Logarithmic" $file_out | tail -1`
          ord_gap=`grep -a -A2 "Ordinary" $file_out | tail -1`
          tot_deg=`grep -a -A1 "degenerate" $file_out | tail -1`
          echo $nspin $V $hz $alpha $T0 $gap_ratio $log_gap $shift_log_gap $ord_gap $tot_deg $iterations >> raw_$output

        done
        echo "" >> raw_$output
        echo ""
      done
    done
  done
  column -t < raw_$output >> temp
  cat > input_sort.txt << *
$nparam
`echo ${cols[@]} | tr ' ' '\n'`
${nums[${cols[$nparam - 1]} - 1]}
*
  ./sort_data.sh temp < input_sort.txt | tee $file_sort
  mv sort_col${cols[$nparam - 1]}_temp $output
  mv temp output/
  
  cat > input_sort.txt << *
$nparam
`echo ${cols[@]:1:$nparam - 1} | tr ' ' '\n'`
${cols[0]}
$block
*
  ./sort_data.sh $output < input_sort.txt | tee $file_sort
done

mv temp $file_out $file_sort raw_$output output/
