
#read file

file=$1

sort -n -k2 -k1 $file -o sort_$file

#Use the following line to delete all blank lines in a given file:
# sed -i '/^$/d' <filename>
sed -i '/^$/d' sort_$file

step=6 #16 #20 #Number of block length plus one ( # of J + 1)
  j=0
for J in $(seq 0.05 0.05 0.30) 0.40 0.50
do
  for V in 0.50
  do
    for hz in 6.00 
    do
      j=$(( $j + 1 ))
      line=$(( $step * $j ))
      #echo $j, $line
      cat > ed_script.ed << *
${line}i

.
w
q
*
      ed sort_$file < ed_script.ed >> ed_out.txt
    done
  done
done

cat sort_$file
