
#read file

file=$1

sort -n -k3 -k1 $file -o sort_$file

#Use the following line to delete all blank lines in a given file:
# sed -i '/^$/d' <filename>
sed -i '/^$/d' sort_$file

step=8 #20 #Number of block length plus one 
  j=0
for param in $(seq 0 0.05 0.30) 0.50 0.79 0.80
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

cat sort_$file
