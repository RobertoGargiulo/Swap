filestring="swap"
gfortran -c mod_exp.f90 mod_genmat.f90 mod_print.f90 $filestring.f90
gfortran mod_exp.o mod_genmat.o mod_print.o $filestring.o -llapack -lblas -fexternal-blas
rm mod_exp.o mod_genmat.o mod_print.o $filestring.o exponentiate.mod genmat.mod printing.mod printing_subrtns.mod
mv a.out $filestring
#./$filestring
gfortran -o mag_avg mag_avg.f90

#: <<'END_COMMENT'
nspin=4
steps=1000
iterations=40
#for nspin in 2 4 6 8
#do
  for kick in 0 #0.00 0.05 0.10 0.40
  do
    for h_coupling in 0.1 #0.00 0.05 0.15 0.30 0.60
    do
      cat > input.txt << *
  $nspin
  $steps
  $iterations
  $h_coupling
  $kick
*
      ./$filestring < input.txt
      ./mag_avg < input.txt
    done
  done
#done
#END_COMMENT
