##!/bin/bash
filestring="test_khemani"

make $filestring

n_threads=1
export OMP_NUM_THREADS=$n_threads 

output="test_Khemani_solvable.txt"
mv $output data/eigen/
file_out="out2.txt"
file_sort="out_sort2.txt"
mv temp $file_out $file_sort raw_$output output/

###Choice of parameters
iterations_2=20480
list_nspin=$(seq 2 2 4)
list_lambda="0.00 1.00" #`echo 9 | awk '{ list=""; val = 0.01; for(i=1;i<=$1;i++) {list = list val " "; val = val * 2}; print list }'`
list_V="1.00"
list_kick="0.00"
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
  #iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  iterations=1
  for kick in $list_kick
  do
    for T0 in $list_T0
    do
      for lambda in $list_lambda
      do
        hx="0.00" #`echo $lambda | awk '{print $1 * 0.1}'`
        hy="0.00" #`echo $lambda | awk '{print $1 * 0.15}'`
        hz=`echo $lambda | awk '{print $1 * 0.45}'`
        for V in $list_V
        do
          cat > input.txt << *
$nspin
$iterations
$T0
$hx
$hy
$hz
$V
$kick
*
  	      ./$filestring < input.txt | tee $file_out
          echo "nspin = $nspin  lambda = $lambda  V = $V  epsilon = $kick  T0 = $T0  hx = $hx hy = $hy hz = $hz"
          echo "iterations = $iterations   n_threads = $n_threads"
          gap_ratio=`grep -a -A2 "Ratio" $file_out | tail -1`
          log_gap=`grep -a -A2 "pi-Logarithmic" $file_out | tail -1`
          shift_log_gap=`grep -a -A2 "Half Shifted Logarithmic" $file_out | tail -1`
          echo $nspin $lambda $V $kick $T0 $hx $hy $hz $gap_ratio $log_gap $shift_log_gap $iterations >> raw_$output
          echo ""
          echo "" >> raw_$output
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
