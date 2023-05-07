##!/bin/bash
filestring="spectrum_Khemani"

make $filestring

n_threads=20
export OMP_NUM_THREADS=$n_threads 

output="Khemani_spectrum_L10_lambda.txt"
mv $output data/eigen/
file_out="out2.txt"
file_sort="out_sort2.txt"
mv temp $file_out $file_sort raw_$output output/

###Choice of parameters
iterations_2=20480
list_nspin=$(seq 4 2 10)
list_lambda="0.01 0.05 0.10 0.50 1.00"
list_V="1.00"
list_kick="0.05"
list_T0="1.00"


nums=()
nums+=(`wc -w <<< $list_nspin`)  
nums+=(`wc -w <<< $list_lambda`)
nums+=(`wc -w <<< $list_V`)
nums+=(`wc -w <<< $list_kick`)
nums+=(`wc -w <<< $list_T0`)

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
      for hx in $list_hx
      do
        for V in $list_V
        do
          for hz in $list_hz
          do
            cat > input.txt << *
$nspin
$iterations
$T0
$hx
$V
$hz
$kick
*
  	        ./$filestring < input.txt | tee $file_out
            echo "nspin = $nspin  hx = $hx  V = $V  hz = $hz  epsilon = $kick  T0 = $T0"
            echo "iterations = $iterations   n_threads = $n_threads"
            gap_ratio=`grep -a -A2 "Ratio" $file_out | tail -1`
            log_gap=`grep -a -A2 "Logarithm" $file_out | tail -1`
            echo $nspin $hx $hy $hz $V $kick $T0 $gap_ratio $log_gap $iterations >> raw_$output
            echo ""
            echo "" >> raw_$output
          done
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
