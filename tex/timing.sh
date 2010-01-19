#!/bin/bash

stop=13
if [ "$1" != "" ]; then
    stop=$1
fi

rm -rf padded unpadded
touch padded unpadded

echo Timing:
for (( i=7; i<=$stop; i++ ))
do
    echo $i
    m=$(asy -c "2^$i")
    echo -e "$m \t $(./a.out $m | grep -A 1 Unpadded | tail -n 1)"| cat >> unpadded
    m=$(asy -c "ceil(2^$i*2/3)")
    echo -e "$m \t $(./a.out $m | grep -A 1 Padded | tail -n 1)" | cat >> padded
done
