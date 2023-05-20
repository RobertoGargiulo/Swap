##!/bin/bash
filestring="test_khemani"

make $filestring

n_threads=10
export OMP_NUM_THREADS=$n_threads 

###Choice of parameters
iterations_2=20480
list_nspin=$(seq 1 1 4)
list_lambda="0.00 1.00" #`echo 9 | awk '{ list=""; val = 0.01; for(i=1;i<=$1;i++) {list = list val " "; val = val * 2}; print list }'`
list_V="1.00"
list_kick="0.00"
list_T0="1.00"

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
        hx=`echo $lambda | awk '{print $1 * 0.1}'`
        hy=`echo $lambda | awk '{print $1 * 0.15}'`
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
          echo ""
        done
      done
    done
  done
done
