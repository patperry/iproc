#!/bin/bash

START=$((1 + 50 * ($1 - 1)))
END=$(($START + 49))

echo "START = ${START}"
echo "END = ${END}"

for i in $(eval echo {$START..$END})
do
    CMD="./tests/enronfit --start=./output/full/full.h5 --boot=$i --resid --output=./output/full/fit$i.h5 > ./output/full/fit$i.out 2> ./output/full/fit$i.err"
    echo $CMD
    eval $CMD
done
