filestring="prova2"
compiler="ifort"
#$compiler -c mod_exp.f90 mod_genmat.f90 mod_print.f90 $filestring.f90 -fopenmp
#$compiler mod_exp.o mod_genmat.o mod_print.o $filestring.o -Wall -Wextra -llapack -lblas -fopenmp
$compiler -c mod_exp.f90 mod_genmat.f90 mod_print.f90 $filestring.f90 -qopenmp
$compiler mod_exp.o mod_genmat.o mod_print.o $filestring.o -llapack -lblas -fopenmp
#$compiler -c mod_exp.f90 mod_genmat.f90 mod_print.f90 $filestring.f90 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
#$compiler mod_exp.o mod_genmat.o mod_print.o $filestring.o -I"${MKLROOT}/include"
rm mod_exp.o mod_genmat.o mod_print.o $filestring.o exponentiate.mod genmat.mod printing.mod printing_subrtns.mod
mv a.out $filestring
#$compiler -o $filestring $filestring.f90 -fopenmp
#./$filestring
#$compiler -o mag_avg mag_avg.f90

#: <<'END_COMMENT'
nspin=6
steps=10  #Time in units of hbar/V
iterations=3
export OMP_NUM_THREADS=20
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
$iterations
$steps
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
