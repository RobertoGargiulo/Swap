file_string="swap"

for j in {1..10000..1}
do
  kill_pid=`ps -u robgarg | grep $file_string | awk '{print $1}' | tr '\n' ' ' && echo`
  #if [ "$kill_pid" != "" ];
  #then
    kill $kill_pid && echo $kill_pid && sleep 0.01
  #fi
done
