
#read file

file=$1

sort -n -k4 -k1 $file -o sort_$file

#Use the following line to delete all blank lines in a given file:
# sed -i '/^$/d' <filename>
sed -i '/^$/d' sort_$file

step=7 #16 #20 #Number of block length plus one ( # of J + 1)
  j=0
for J in 0.10 #$(seq 0.05 0.05 0.30) #$(seq 0.01 0.01 0.10; seq 0.20 0.10 1.00) #0.10
do
  for V in 0.50
  do
    for hz in 0.01 2.00 $(seq 4.00 4.00 16.00) #6.00 # 0.01 $( seq 0.20 0.20 2.00; seq 4.00 4.00 16.00)
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
