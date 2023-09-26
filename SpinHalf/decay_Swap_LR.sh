##!/bin/bash
filestring="decay_swap_LR"

make $filestring

n_threads=10
export OMP_NUM_THREADS=$n_threads 

output="Swap_LR_decay_kick0_L12.txt"
sorting=false
file_out="out.txt"
file_sort="out_sort.txt"
if [ "$sorting" = true ] ; then
  mv $output data/eigen/
  mv temp $file_out $file_sort raw_$output output/
fi

###Choice of parameters
n_disorder_2=20480
steps="1000000"
list_nspin=$(seq 10 2 12)
list_J="0.15 0.08 0.065 0.01" #"2.0 0.5 0.1 0.05" #0.001 0.0001 0.00001"
list_V="3.00"
list_hz="16.00"
list_alpha="0.50 3.00" #1.00 10.00
list_T0="1.00"
list_kick="0.00"


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
  #n_disorder=20
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
              if (( $(echo "$J >= 0.25" |bc -l) )); then
                  n_pow_periods=0
              #elif (( $(echo "$J >= 0.1 && $J < 0.25" |bc -l) )); then
              #  if [ $nspin -le 4 ]; then
              #    n_pow_periods=$((nspin/2 - 1))
              #  elif [ $nspin -gt 4 ] && [ $nspin -le 10 ]; then
              #    n_pow_periods=$((2*nspin - 4))
              #  elif [ $nspin -gt 10 ]; then
              #    n_pow_periods=$((nspin))
              #    fi
              #elif (( $(echo "$J < 0.1" |bc -l) )); then
              #  if [ $nspin -le 4 ]; then
              #    n_pow_periods=$((nspin+1))
              #  elif [ $nspin -gt 4 ] && [ $nspin -le 10 ]; then
              #    n_pow_periods=$((nspin+2))
              #  elif [ $nspin -gt 10 ]; then
              #    n_pow_periods=$((nspin+2))
              #  fi
              elif (( $(echo "$J > 0.05 && $J < 0.25" |bc -l) )); then
                n_pow_periods=$(awk -v J="${J}" -v L="${nspin}" -v steps="${steps}" 'BEGIN { {res=int(log( 0.42 * (3.1*J)**(-0.38*L) )/log(2))+1 } { printf "%d", (res>=0)?res:0} }');
              elif (( $(echo "$J > 0.01 && $J <= 0.05" |bc -l) )); then
                n_pow_periods=$(awk -v J="${J}" -v L="${nspin}" -v steps="${steps}" 'BEGIN { {res=int(log( 0.42 * (3.1*J)**(-0.38*L) )/log(2))+2 } { printf "%d", (res>=0)?res:0} }');
                if [ $nspin -ge 8 ]; then
                  n_pow_periods=$((n_pow_periods+1))
                fi
              elif (( $(echo "$J <= 0.01" |bc -l) )); then
                n_pow_periods=$(awk -v J="${J}" -v L="${nspin}" -v steps="${steps}" 'BEGIN { {res=int(log( 0.42 * (3.1*J)**(-0.38*L) )/log(2))+3 } { printf "%d", (res>=0)?res:0} }');
                if [ $nspin -ge 8 ]; then
                  n_pow_periods=$((n_pow_periods+1))
                fi
              fi
              n_periods=$(( 2**(n_pow_periods) ))
              total_time=$(( n_periods * steps ))
              cat > input.txt << *
$nspin
$n_disorder
$steps
$n_pow_periods
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
              echo "n_disorder = $n_disorder steps = $steps n_pow_periods = $n_pow_periods n_periods = $n_periods total_time = $total_time n_threads = $n_threads"
              echo ""
              echo "Executing file: "
              ./$filestring < input.txt | tee $file_out
              #echo "Recap and Output: "
              #echo "nspin = $nspin  J = $J  V = $V  hz = $hz  epsilon = $kick  alpha = $alpha  T0 = $T0"
              #echo "n_disorder = $n_disorder  steps = $steps n_pow_periods = $n_pow_periods total_time = $total_time n_threads = $n_threads"
              decay=`grep -a -A2 "Decay Time" $file_out | tail -1`
              echo $nspin $J $V $hz $kick $alpha $T0 $decay $n_pow_periods $n_periods $total_time $steps $n_disorder >> raw_$output
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
