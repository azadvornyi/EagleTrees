#!/bin/bash

python3 noh.py 100 > holder

while read here; do

    for (( c=0; c<${here}; c++ )) do
          python3 sat_tracker.py $c 100 #insert [sat_tracker.py] instead of dummy1.py
      done

done < holder




while read name; do

    for i in "${name[@]}"; do
        A="$(cut -d' ' -f1 <<<"$i")"
        B="$(cut -d' ' -f2 <<<"$i")"
        for (( k=0; k < $B; k++ )) do
          python3 simpair.py $A $k 100
        done
    done

done < hostid_satnumber_100.txt



SUBJECT="Calculations are finished"
TO="a.zadvornyi@student.rug.nl"
MESSAGE="message.txt"

/usr/bin/mail -s "$SUBJECT" "$TO" < $MESSAGE