##!/bin/bash
filestring="spectrum_swap_LR"

make $filestring

n_threads=10
export OMP_NUM_THREADS=$n_threads 

output="Swap_LR_spectrum_L12_J_alpha_large_N_dis.txt"
mv $output data/eigen/
file_out="out.txt"
file_sort="out_sort.txt"
mv temp $file_out $file_sort raw_$output output/

###Choice of parameters
iterations_2=20480
list_nspin=$(seq 4 2 12)
list_J="0.001 0.01 0.1"
list_V="3.00"
list_hz="16.00"
list_kick="0.001 0.01 0.05"
list_alpha="0.50 1.00 3.00 10.00" #20.00
list_T0="1.00"

n1=`wc -w <<< $list_nspin`
n2=`wc -w <<< $list_J`
n3=`wc -w <<< $list_V`
n4=`wc -w <<< $list_hz`
n5=`wc -w <<< $list_kick`
n6=`wc -w <<< $list_alpha`
n7=`wc -w <<< $list_T0`

list_nparam="$n1 $n2 $n3 $n4 $n5 $n6 $n7"
echo "list_param = "$list_nparam
nums=($list_nparam)
echo "nums = ${nums[@]}"
cols=`echo "${nums[@]}" | awk '{ indexes=""; for(i=1;i<=NF;i++) {if($i>1) {indexes=indexes i " "}} print indexes }'`
echo "cols = "$cols
cols=($cols)
nparam=${#cols[@]}


block=0
for nspin in $list_nspin
do
  block=$(( $block + 1 ))
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
          #for idx in {0..10..1}
          for hz in $list_hz
          do
            #hz=`echo $(printf "%.2f" $(echo "0.04 * 2^($idx)" | bc) )`
            for alpha in $list_alpha
            do
              cat > input.txt << *
$nspin
$iterations
$T0
$J
$V
$hz
$kick
$alpha
*
              file_out="out.txt"
  	          ./$filestring < input.txt | tee $file_out
              echo "nspin = $nspin  J = $J  V = $V  hz = $hz  epsilon = $kick  alpha = $alpha  T0 = $T0"
              echo "iterations = $iterations   n_threads = $n_threads"
              gap_ratio=`grep -a -A2 "Ratio" $file_out | tail -1`
              log_gap=`grep -a -A2 "Logarithm" $file_out | tail -1`
              echo $nspin $J $V $hz $kick $alpha $T0 $gap_ratio $log_gap $iterations >> raw_$output
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
  cat > input_sort.txt << *
$nparam
`echo ${cols[@]} | tr ' ' '\n'`
${nums[${cols[$nparam - 1]} - 1]}
*
  #echo "input_sort.txt: "
  #cat input_sort.txt
  #echo "cols = "${cols[@]}
  #echo "nums = "${nums[@]}
  #echo "nparam = "$nparam
  #echo "nums[cols[nparam]] = "${nums[${cols[$nparam - 1]} - 1]}
  #./sort_data.sh temp < input_sort.txt | tee $file_sort
  #mv sort_col${cols[$nparam - 1]}_temp $output
  #mv temp output/

  cat > input_sort.txt << *
$nparam
`echo ${cols[@]:1:$nparam - 1} | tr ' ' '\n'`
${cols[0]}
$block
*
  #echo "input_sort.txt: "
  #cat input_sort.txt
  #echo "cols[2] = "${cols[1]}
  #echo "cols[2:nparam] = ${cols[@]:1:$nparam - 1}"
  #echo "cols[1] = "${cols[0]}
  
  ./sort_data.sh $output < input_sort.txt | tee $file_sort
done

mv temp $file_out $file_sort raw_$output output/
