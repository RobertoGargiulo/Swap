

#Use the following line to delete all blank lines in a given file:
# sed -i '/^$/d' <filename>


file="sort_Swap_decay_times.txt"

  j=0
  for J in $(seq 0.02 0.02 0.10)
  do
      for V in 0.25 0.50
      do
        for hz in 3.00 6.00 12.00
        do
        j=`echo $j | awk '{print $1+1}'`
        #echo $j
        line=`echo $j | awk '{print 6*$1}'` 
        cat > ed_script.ed << *
${line}i

.
w
q
*
        ed $file < ed_script.ed >> out_ed.txt
      done
    done
  done
