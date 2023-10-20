file_string=$1

for j in {1..10000..1}
do
  kill_pid=`ps -u gargiulo | grep $file_string | awk '{print $1}' | tr '\n' ' ' && echo`
  #if [ "$kill_pid" != "" ];
  #then
    kill -KILL $kill_pid && echo $kill_pid && sleep 0.3
  #fi
done
