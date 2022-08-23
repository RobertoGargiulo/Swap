##!/bin/bash
filestring="swap_decay"

make $filestring

n_threads=30
#total_time=1000000
export OMP_NUM_THREADS=$n_threads 
output="Swap_decay_times_varying_hz.txt"
#mv $output ../Trash

iterations_2=5120
for nspin in 12 #8 10 #12
do
  iterations=`echo $iterations_2 $nspin | awk '{print 2**(-$2/2+1)*$1}'`
  for kick in 0.00 #$(seq 0.00 0.05 0.20) 
  do
    for period in 1.00 
    do
      for J in 0.10 # $(seq 0.05 0.05 0.30) 
      do
        for V in 0.25 #$(seq 0.00 0.40 1.60) 
        do
          for hz in 12.00 16.00 #0.01 2.00 $(seq 4.00 4.00 16.00)
          do
            if (( $(echo "$hz <= 4.00" |bc -l) )); then
              if [ $nspin -le 4 ]; then
                total_time=$((10**4))
                n_periods=1
              elif [ $nspin -gt 4 ] && [ $nspin -le 8 ]; then
                total_time=$((10**5))
                n_periods=1
              elif [ $nspin -gt 8 ] && [ $nspin -le 12 ]; then
                total_time=$((5*10**5))
                n_periods=10
              fi
            elif (( $(echo "$hz > 4.00 && $hz <= 8.00" |bc -l) )); then
              if [ $nspin -le 4 ]; then
                total_time=$((10**(nspin/2+3)))
                n_periods=1
              elif [ $nspin -gt 4 ] && [ $nspin -le 8 ]; then
                total_time=$((10**(nspin/2+3)))
                n_periods=10
              elif [ $nspin -gt 8 ] && [ $nspin -le 12 ]; then
                total_time=$((10**(nspin/2+3)/2))
                n_periods=500
              fi
            elif (( $(echo "$hz > 8.00" |bc -l) )); then
              if [ $nspin -le 4 ]; then
                total_time=$((10**(nspin/2+3)))
                n_periods=1
              elif [ $nspin -gt 4 ] && [ $nspin -le 8 ]; then
                total_time=$((10**(nspin/2+4)))
                n_periods=100
              elif [ $nspin -gt 8 ] && [ $nspin -le 12 ]; then
                total_time=$((5*10**(nspin/2+5)))
                n_periods=5000000
              fi
            fi
            steps=`echo $total_time $period $n_periods | awk '{print $1/$2/$3}'`

            cat > input.txt << *
$nspin
$iterations
$steps
$n_periods
$period
$J
$V
$hz
$kick
*
            file_out="out5.txt"
  	        ./$filestring < input.txt | tee $file_out
            echo "nspin = $nspin "
            echo "J = $J  V = $V  hz = $hz  epsilon = $kick"
            echo "period = $period  n_periods = $n_periods  steps = $steps   total_time = $total_time "
            echo "iterations = $iterations   n_threads = $n_threads "
            decay_times=`grep -a -A1 "Decay Time" $file_out | tail -1`
            echo $nspin $J $V $hz $kick $period $decay_times $iterations $total_time $n_periods $steps >> $output
          done
          echo "" >> $output
          echo ""
        done
      done
    done
  done
done
