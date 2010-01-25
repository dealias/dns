#!/bin/bash

start=2
stop=6
if [ "$1" != "" ]; then
    stop=$1
fi

rm -rf padded unpadded
touch padded unpadded

echo Timing:
for (( i=$start; i<=$stop; i++ ))
do
    echo $i
    m=$(asy -c "2^$i")
    echo -e "$m \t $(./a.out $m 0 | grep -A 1 Unpadded | tail -n 1)"| cat >> unpadded
    m=$(asy -c "2^$i")
    echo -e "$m \t $(./a.out $m 1 | grep -A 1 Padded | tail -n 1)" | cat >> padded
done
