#!/bin/bash

rm -rf padded unpadded
touch padded unpadded
echo Timing:
for i in {6..14}
do
    echo $i
    n=$(asy -c "2^$i")
    echo -e "$n \t $(./a.out $n | grep -A 1 Unpadded | tail -n 1)"| cat >> unpadded
    n=$(asy -c "ceil(2^$i*2/3)")
    echo -e "$n \t $(./a.out $n | grep -A 1 Padded | tail -n 1)" | cat >> padded
done
