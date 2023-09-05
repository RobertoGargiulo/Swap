
#read file

file=$1

#print to output for comparison with sorted file
#echo "Original file: "
#cat $file
#echo ""

echo "Number of system sizes L (aka block size): "
read bsize

echo "Number of points (J, epsilon, ...) BUT L (an upper bound also works): "
read npoints

output="sort_w_L_${file}"
mv $output output/

step=$((bsize + 1)) #Step size is equal to the block size plus the one blank line


#sort files numerically according first to column $1 and column $2
#sort -n $ordering $file -o $output 
sort -n -k6 -k2 -k3 -k4 -k5 -k7 -k1 $file -o $output


#Delete all blank lines in a given file
sed -i '/^$/d' $output

#Re-insert them at the appropriate lines
j=0
line=0
nlines=`wc -l < $file`
while : ; do
  j=$(( $j + 1 ))
  line=$(( $step * $j ))

  cat > ed_script.ed << *
${line}i

.
w
q
*
  #Insert blank line at position $step * $j using 'ed'
  ed $output < ed_script.ed >> ed_out.txt
  nlines=`wc -l < $output`
  #echo $j, $line, $nlines
  
  #Condition for while-loop end
  [[ $line -lt $nlines ]] || break
done

#print to output to see result
#echo "Re-ordered file: "
#cat $output
