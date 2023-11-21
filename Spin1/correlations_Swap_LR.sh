##!/bin/bash
filestring="correlations_swap_LR"

make $filestring

ulimit -s unlimited
n_threads=50
stack=6G
export OMP_NUM_THREADS=$n_threads 
export OMP_STACKSIZE=$stack
sleep_time=0.0

output="Swap_LR_correlations_L10_large_J.txt"
sorting=true
file_out="out_3.txt"
file_sort="out_sort_3.txt"
if [ "$sorting" = true ] ; then
  mv $output data/eigen/
  mv temp $file_out $file_sort raw_$output output/
fi

###Choice of parameters
n_disorder_2=5120 #20480
list_nspin="10" #$(seq 4 2 10)
list_J=`echo $(seq 8 1 10) | tr ' ' '\n' | awk '{print 0.005*2**($1)}'`
list_V="3.00"
list_hz="16.00"
list_kick="0.0"
list_alpha="0.50 3.00"
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
  n_disorder=`echo $n_disorder_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  #n_disorder=10
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
$n_disorder
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
              echo "n_disorder = $n_disorder n_threads = $n_threads omp_stack_size = $stack"
              echo ""
              echo "Executing file: "
              sleep $sleep_time 
  	          ./$filestring < input.txt | tee $file_out
              corr=`grep -a -A1 "Total average correlations:" $file_out | tail -1`
              norm_corr=`grep -a -A1 "Total average normalized correlations:" $file_out | tail -1`
              echo $nspin $J $V $hz $kick $alpha $T0 $corr $norm_corr $n_disorder >> raw_$output
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
