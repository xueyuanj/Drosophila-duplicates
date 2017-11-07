#!/bin/bash

c=1; d=1; mkdir -p dir_${d}

for jpg_file in *.fasta
do
        if [ $c -eq 50 ]
        then
                d=$(( d + 1 )); c=0; mkdir -p dir_${d}
        fi
        mv "$jpg_file" dir_${d}/
        c=$(( c + 1 ))
done