#!/bin/bash




python3 lookingforbugs.py 25


while read here; do
    for i in "${here[@]}"; do
          A="$(cut -d' ' -f1 <<<"$i")"
          B="$(cut -d' ' -f2 <<<"$i")"
          C="$(cut -d' ' -f3 <<<"$i")"
          D="$(cut -d' ' -f4 <<<"$i")"
          E="$(cut -d' ' -f5 <<<"$i")"
          F="$(cut -d' ' -f6 <<<"$i")"
          python3 sat_tracker.py $A $B $C $D $E $F 25 #insert [sat_tracker.py] instead of dummy1.py
      done
done < hostid_rvir_33_25.txt

while read name; do
    for i in "${name[@]}"; do
        A="$(cut -d' ' -f1 <<<"$i")"
        B="$(cut -d' ' -f2 <<<"$i")"
        C="$(cut -d' ' -f3 <<<"$i")"
        D="$(cut -d' ' -f4 <<<"$i")"
        E="$(cut -d' ' -f5 <<<"$i")"
        F="$(cut -d' ' -f6 <<<"$i")"
        G="$(cut -d' ' -f7 <<<"$i")"
        for (( k=0; k < $G; k++ )) do

            python3 simpair.py $A $B $C $D $E $F $k 25

        done
    done

done < hostid_satnumber_rvir_33_25.txt

echo "done..."

SUBJECT="Calculations are finished"
TO="a.zadvornyi@student.rug.nl"
MESSAGE="message.txt"

/usr/bin/mail -s "$SUBJECT" "$TO" < $MESSAGE