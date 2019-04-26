#!/bin/bash



for (( c=0; c<36; c++ ))
do
   python3 dummy1.py $c
done

host=0
while read name; do

    for (( k=0; k<=$name; k++ ))
    do
        python3 simpair.py $host $k
    done
    host=$(($host+1))

done < satnumber.txt