filestring="prova"
compiler="gfortran"
export OMP_NUM_THREADS=20
$compiler -c mod_exp.f90 mod_genmat.f90 mod_print.f90 $filestring.f90
$compiler mod_exp.o mod_genmat.o mod_print.o $filestring.o -Wall -Wextra -llapack -lblas -fexternal-blas -fopenmp
rm mod_exp.o mod_genmat.o mod_print.o $filestring.o exponentiate.mod genmat.mod printing.mod printing_subrtns.mod
mv a.out $filestring
#./$filestring
$compiler -o mag_avg mag_avg.f90

#: <<'END_COMMENT'
nspin=8
steps=10  #Time in units of hbar/V
iterations=4
time_step=0.05
#for nspin in 2 4 6 8
#do
  for kick in 0.00 #0.05 0.10 0.20 0.40
  do
    for J in 1 #0 0.5 1 1.5 2
    do
      for V in 2 #0 0.5 1 1.5 2 3
      do
        for hx in 0 #0.15 0.30 0.60 1.00
        do
          for hz in 3 #0 1 2.5 3 4 6
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
            #./mag_avg < input.txt
          done
        done
      done
    done
  done
#done
#END_COMMENT
