#!/usr/bin/env bash

python3 noh.py > holder



while read here; do

    for (( c=0; c<${here}; c++ )) do
          python3 dummy1.py $c  #insert [sat_tracker.py] instead of dummy1.py
      done

done < holder

