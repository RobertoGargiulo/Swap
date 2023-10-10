##!/bin/bash
filestring="dynamics_swap_LR"

ulimit -s unlimited
make $filestring

n_threads=32
stack=6G
export OMP_NUM_THREADS=$n_threads 
export OMP_STACKSIZE=$stack

output="Swap_LR_dynamics_long.txt"
file_out="out.txt"
mv $file_out output/

###Choice of parameters
n_disorder_2=10240
steps="1000000"
list_nspin=$(seq 4 2 12)
list_J="0.1 0.01"
list_V="3.00"
list_hz="16.00"
list_kick="0.00"
list_alpha="0.50 3.00" #1.00 10.00
list_T0="1.00"

for nspin in $list_nspin
do
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
              if (( $(echo "$J >= 0.25" |bc -l) )); then
                  n_pow_periods=0
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
            done
            echo ""
          done
        done
      done
    done
  done
done
