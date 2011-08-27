#!/bin/bash

for i in {1..10}
do
    CMD="nohup ./boot.sh $i > run$i.out 2> run$i.err &"
    echo $CMD
    eval $CMD
done
