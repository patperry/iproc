#!/bin/bash

START=$((1 + 50 * ($1 - 1)))
END=$(($START + 49))

echo "START = ${START}"
echo "END = ${END}"

for i in $(eval echo {$START..$END})
do
    CMD="./tests/enronfit --start=output/fit-dynamic.json --boot=$i > ./output/boot-dynamic/fit$i.json 2> ./output/boot-dynamic/fit$i.err"
    echo $CMD
    eval $CMD
done
