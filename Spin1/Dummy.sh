filestring="dummy"
make $filestring
n_threads=20
export OMP_NUM_THREADS=$n_threads 

for i in {1..10}
do
  echo "i = ", $i
  ./$filestring
  echo ""
  echo ""
done
