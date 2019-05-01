#!/usr/bin/env bash

start=`date +%s`



#python3 1h+1s.py 1 1 & python3 1h+1s.py 1 2  & python3 1h+1s.py 1 3  & python3 1h+1s.py 1 4  && fg
#parallel ::: 1h+1s.py 1 1  python3 1h+1s.py 1 2   python3 1h+1s.py 1 3   python3 1h+1s.py 1 4


python3 1h+1s.py 1 1 &
P1=$!
python3 1h+1s.py 1 2 &
P2=$!
python3 1h+1s.py 1 3 &
P3=$!
#python3 1h+1s.py 1 4 &
#P4=$!
#python3 1h+1s.py 1 5 &
#P5=$!
#python3 1h+1s.py 1 6 &
#P6=$!
wait $P1 $P2 $P3 #$P4 #$P5 #$P6

end=`date +%s`

runtime=$((end-start))

echo "$runtime s takes to run 5 scripts in parallel-----------------------------"