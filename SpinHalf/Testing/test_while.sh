j=0
if [ "$j" -gt 10 ]; then
  echo $j
else
  echo "If statement is working"
fi
while : ; do
    j=$(( ${j} + 1 ))
    echo "Loop" ${j}
    [ "$j" -lt 10 ] || break
done

