filestring="prova3"

make $filestring

#: <<'END_COMMENT'
read nspin
steps=20
iterations=1
time_step=0.005
#for nspin in 2 4 6 8
#do
  for kick in 0.00 #0.05 0.10 0.20 0.40
  do
    for J in 2 #1 1.5 2
    do
      for V in 1 #0 0.5 1 1.5 2 3
      do
        for hz in 1 #0 1 2.5 3 4 6
        do
          cat > input.txt << *
$nspin
$iterations
$steps
$time_step
$J
$V
$hz
$kick
*
          ./$filestring < input.txt
          #./mag_avg < input.txt
        done
      done
    done
  done
#done
#END_COMMENT
