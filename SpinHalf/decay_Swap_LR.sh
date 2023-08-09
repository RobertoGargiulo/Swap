##!/bin/bash
filestring="decay_kick0_swap_LR"

make $filestring

n_threads=8
export OMP_NUM_THREADS=$n_threads 

output="Swap_LR_decay_kick0_L12.txt"
sorting=true
file_out="out.txt"
file_sort="out_sort.txt"
if [ "$sorting" = true ] ; then
  mv $output data/eigen/
  mv temp $file_out $file_sort raw_$output output/
fi

###Choice of parameters
n_disorder_2=20480
steps="1000000"
list_nspin=$(seq 4 2 12)
list_J="1.0 0.1 0.01 0.001 0.0001 0.00001"
list_V="3.00"
list_hz="16.00"
list_alpha="0.50 3.00" #1.00 10.00
list_T0="1.00"


nums=()
nums+=(`wc -w <<< $list_nspin`)
nums+=(`wc -w <<< $list_J`)
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
  n_disorder=`echo $n_disorder_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  n_disorder=20
  for T0 in $list_T0
  do
    for J in $list_J
    do
      if (( $(echo "$J > 0.01" |bc -l) )); then
        if [ $nspin -le 4 ]; then
          total_time=$((10**4))
          n_periods=1
        elif [ $nspin -gt 4 ] && [ $nspin -le 8 ]; then
          total_time=$((10**5))
          n_periods=1
        elif [ $nspin -gt 8 ] && [ $nspin -le 12 ]; then
          total_time=$((5*10**5))
          n_periods=11
        fi
      elif (( $(echo "$J > 0.001 && $J <= 0.01" |bc -l) )); then
        if [ $nspin -le 4 ]; then
          total_time=$((10**(nspin/2+3)))
          n_periods=1
        elif [ $nspin -gt 4 ] && [ $nspin -le 8 ]; then
          total_time=$((10**(nspin/2+3)))
          n_periods=11
        elif [ $nspin -gt 8 ] && [ $nspin -le 12 ]; then
          total_time=$((10**(nspin/2+3)/2))
          n_periods=101
        fi
      elif (( $(echo "$J < 0.001" |bc -l) )); then
        if [ $nspin -le 4 ]; then
          total_time=$((10**(nspin/2+3)))
          n_periods=1
        elif [ $nspin -gt 4 ] && [ $nspin -le 8 ]; then
          total_time=$((10**(nspin/2+4)))
          n_periods=101
        elif [ $nspin -gt 8 ] && [ $nspin -le 12 ]; then
          total_time=$((5*10**(nspin/2+5)))
          n_periods=501
        fi
      fi
      steps=`echo $total_time $T0 $n_periods | awk '{print int($1/$2/$3)}'`

      for V in $list_V
      do
        for hz in $list_hz
        do
          for alpha in $list_alpha
          do
            cat > input.txt << *
$nspin
$n_disorder
$n_periods
$steps
$T0
$J
$V
$hz
$alpha
*
            ./$filestring < input.txt | tee $file_out
            echo "nspin = $nspin  J = $J  V = $V  hz = $hz  epsilon = $kick  alpha = $alpha  T0 = $T0"
            echo "n_disorder = $n_disorder  steps = $steps n_periods = $n_periods total_time = $total_time n_threads = $n_threads"
            decay=`grep -a -A2 "Decay Time" $file_out | tail -1`
            echo $nspin $J $V $hz $kick $alpha $T0 $decay $n_periods $total_time $steps $n_disorder >> raw_$output
          done
          echo ""
          echo "" >> raw_$output
        done
        #echo "" >> $output
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
