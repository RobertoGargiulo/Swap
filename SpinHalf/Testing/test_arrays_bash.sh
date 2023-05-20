#function my_function() {
#  local -a nums=("$@")
#  local -a indexes=($(echo "${!nums[@]}" | awk '{for(i=1;i<=NF;i++)if($i>=1) j=0; j++; a[j]=i} END {for(i in a)print i-1}'))
#  echo "Sum of positive integers: $(echo "${nums[@]}" | tr ' ' '\n' | awk '$1>0{sum+=$1}END{print sum}')"
#  echo "Indexes of positive integers: ${indexes[@]}"
#}

read list
#echo $list
#my_function $list

#function sum_and_index_greater_than_zero() {
#  local -a nums=("$@") indexes=()
#  for i in "${!nums[@]}"; do ((nums[i]>0)) && indexes+=("$i") && ((sum+=nums[i])); done
#  echo "Sum of positive integers: $sum"
#  echo "Indexes of positive integers: ${indexes[@]}"
#}


#function sum_and_index_greater_than_zero() {
#  local -a nums=("$@")
#  echo "${nums[@]}" | awk '{sum=0; indexes=""; for(i=1;i<=NF;i++) {if($i>0) {sum+=$i; indexes=indexes i-1 " "}} print "Sum of positive integers: " sum; print "Indexes of positive integers: " indexes}'
#}
#
#sum_and_index_greater_than_zero $list
#
echo "List = " ${list}
#indx=()
nums=("$list")
echo $nums
echo ${nums[0]}
indx=`echo "${nums[@]}" | awk '{ indexes=""; for(i=1;i<=NF;i++) {if($i>1) {indexes=indexes i " "}} print indexes }'`
#nparam=`echo "${nums[@]}" | awk '{ sum=0; for(i=1;i<=NF;i++) {if($i>1) {sum+=1}} print sum }'`
indx=($indx)  #`echo $indx | tr ' ' '\n'`
nparam=${#indx[@]}
#echo "${nums[@]}" | awk '{sum=0; indexes=""; for(i=1;i<=NF;i++) {if($i>0) {sum+=$i; indexes=indexes i-1 " "}} print "Sum of positive integers: " sum; print "Indexes of positive integers: " indexes}'

echo "indx $indx"
echo "n_param $nparam"

nums=($nums)

cat > input.txt << *
$nparam
`echo ${indx[@]} | tr ' ' '\n'`
${nums[${indx[$nparam - 1]} - 1]}
*
cat input.txt

echo "End file"

cat > input.txt << *
$nparam
`echo ${indx[@]:1:$nparam - 1} | tr ' ' '\n'`
${indx[0]}
${nums[${indx[0]} - 1]}
*
cat input.txt


echo ${nums[@]}
for i in "${!nums[@]}"; do nums[i]=$((${nums[i]}*2)); done
echo ${nums[@]}



