##!/bin/bash
filestring="spectrum_swap_LR"

make $filestring

ulimit -s unlimited
n_threads=40
stack=6G
export OMP_NUM_THREADS=$n_threads 
export OMP_STACKSIZE=$stack

output="Swap_LR_spectrum_non_zero_J_kick.txt"
sorting=true
file_out="out2.txt"
file_sort="out_sort.txt"
if [ "$sorting" = true ] ; then
  mv $output data/eigen/
  mv temp $file_out $file_sort raw_$output output/
fi

###Choice of parameters
iterations_2=40960
list_nspin=$(seq 4 2 12)
list_J="0.01"
list_V="3.00"
list_hz="16.00"
list_kick="0.01 0.005 0.001 0.0005" #0.001 0.01 0.10 0.50
list_alpha="0.50 3.00" #1.00 10.00
list_T0="1.00"


nums=()
nums+=(`wc -w <<< $list_nspin`)
nums+=(`wc -w <<< $list_J`)
nums+=(`wc -w <<< $list_V`)
nums+=(`wc -w <<< $list_hz`)
nums+=(`wc -w <<< $list_kick`)
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
  #####################iterations=`echo $iterations_2 $nspin | awk '{print int( 2**(-$2/2+1)*$1 )}'`##########
  iterations=`echo $iterations_2 $nspin $n_threads | awk '{print int( 1.6**(-$2/2+1)*$1/$3 ) * $3}'`
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
$T0
$J
$V
$hz
$kick
$alpha
*
  	          ./$filestring < input.txt | tee $file_out
              sleep 0.5
              echo "nspin = $nspin  J = $J  V = $V  hz = $hz  epsilon = $kick  alpha = $alpha  T0 = $T0"
              echo "iterations = $iterations   n_threads = $n_threads"
              gap_ratio=`grep -a -A2 "Gap Ratio" $file_out | tail -1`
              log_gap=`grep -a -A2 "pi-Logarithmic" $file_out | tail -1`
              shift_log_gap=`grep -a -A2 "Half Shifted Logarithmic" $file_out | tail -1`
              echo $nspin $J $V $hz $kick $alpha $T0 $gap_ratio $log_gap $shift_log_gap $iterations >> raw_$output
            done
            echo ""
            echo "" >> raw_$output
          done
          #echo "" >> $output
        done
      done
    done
  done
  if [ "$sorting" = true ] ; then
  column -t < raw_$output >> temp

    cat > input_sort.txt << *
$nparam
`echo ${cols[@]:1:$nparam} | tr ' ' '\n'`
${cols[0]}
$block
*
    ./sort_data.sh temp < input_sort.txt | tee $file_sort
    mv sort_col${cols[0]}_temp $output
    mv temp output/

    for i in $(seq 2 1 $(( $nparam - 1)) )
    do
      cat > input_sort.txt << *
$nparam
`echo ${cols[@]:$i:$nparam} | tr ' ' '\n'`
`echo ${cols[@]:0:$i} | tr ' ' '\n'`
${nums[${cols[$i - 1]} - 1]}
*
      ./sort_data.sh $output < input_sort.txt | tee $file_sort
    done

    i=$nparam
    cat > input_sort.txt << *
$nparam
`echo ${cols[@]:0:$i} | tr ' ' '\n'`
${nums[${cols[$i - 1]} - 1]}
*
    ./sort_data.sh $output < input_sort.txt | tee $file_sort
  
  fi

done

if [ "$sorting" = true ] ; then
  mv temp $file_out $file_sort raw_$output output/
fi
