for j in $(seq 1 1 $1)
do
  k=$(( $1 - $j + 1  ))
  #ssh robgarg@qmat.fisica.unina.it "ls -t /home/robgarg/Tesi/Swap/data/eigen | head -$k | tail -1"
  #scp robgarg@qmat.fisica.unina.it:$(ssh robgarg@qmat.fisica.unina.it "ls -t /home/robgarg/Tesi/Swap/data/eigen/*$2* | head -$k | tail -1") data/eigen/
  server="robgarg@qmat.fisica.unina.it"
  dir_remote="/home/robgarg/Tesi/Swap/data/eigen/*$2*"
  dir_local="data/eigen/"
  scp $server:$(ssh $server "ls -t $dir_remote | head -$k | tail -1") $dir_local
done
#/home/robgarg/Tesi/Swap/data/eigen/
