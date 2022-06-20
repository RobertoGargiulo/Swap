filestring="swap"
gfortran -c mod_exp.f90 mod_genmat.f90 mod_print.f90 $filestring.f90
gfortran mod_exp.o mod_genmat.o mod_print.o $filestring.o -llapack -lblas -fexternal-blas
rm mod_exp.o mod_genmat.o mod_print.o $filestring.o exponentiate.mod genmat.mod printing.mod printing_subrtns.mod
mv a.out $filestring
#./$filestring
gfortran -o mag_avg mag_avg.f90

#: <<'END_COMMENT'
nspin=6
steps=20000  #Time in units of hbar/V
iterations=40
time_step=0.2
#for nspin in 2 4 6 8
#do
  for kick in 0.00 #0.05 0.10 0.20 0.40
  do
    for J in 0 0.5 1 2 3
    do
      for V in 1 #0 0.5 1 1.5 2 3
      do
        for hx in 0 #0.15 0.30 0.60 1.00
        do
          for hz in 0 0.5 1 2 3 4 6
          do
            cat > input.txt << *
$nspin
$steps
$iterations
$time_step
$J
$V
$hx
$hz
$kick
*
            ./$filestring < input.txt
            ./mag_avg < input.txt
          done
        done
      done
    done
  done
#done
#END_COMMENT
