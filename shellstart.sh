#!/bin/bash

python3 noh.py > holder

while read here; do

    for (( c=0; c<${here}; c++ )) do
          python3 dummy1.py $c  #insert [sat_tracker.py] instead of dummy1.py
      done

done < holder




while read name; do

    for i in "${name[@]}"; do
        A="$(cut -d' ' -f1 <<<"$i")"
        B="$(cut -d' ' -f2 <<<"$i")"


        for (( k=0; k < $B; k++ )) do
          python3 simpair.py $A $k
        done
    done

done < hostid_satnumber.txt