#!/bin/bash

START=$((1 + 50 * ($1 - 1)))
END=$(($START + 49))

echo "START = ${START}"
echo "END = ${END}"

for i in $(eval echo {$START..$END})
do
    CMD="./tests/enronfit --start=output.h5 --boot=$i --output=./output/boot/fit$i.h5 > ./output/boot/fit$i.out 2> ./output/boot/fit$i.err"
    echo $CMD
    eval $CMD
done
