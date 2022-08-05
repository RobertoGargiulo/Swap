iterations=80
n_threads=1


#Use the following line to delete all blank lines in a given file:
# sed -i '/^$/d' <filename>


export OMP_NUM_THREADS=$n_threads 
for nspin in 12 #6 8 10 12 14 #16 #{2..18..2}
do
  j=0
  for J in 0.50
  do
      for V in $( seq 0.00 0.05 3.00 )
      do
        j=`echo $j | awk '{print $1+1}'`
        #echo $j
        line=`echo $j | awk '{print 122*$1}'` 
        cat > ed_script.ed << *
${line}i

.
w
q
*
        file="data/phases/PT_Sz0_DENSE_MBL_hz_Disorder_AVG_Gap_Ratio_nspin$nspin\_iterations$iterations\_J$J\_up_to_V3_hz6.txt"
        ed $file < ed_script.ed >> out_ed.txt
      done
    done
  done
