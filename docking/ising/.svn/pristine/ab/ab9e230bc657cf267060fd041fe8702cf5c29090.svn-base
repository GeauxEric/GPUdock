#!/bin/bash

mydir=output_`date "+%Y%m%d_%H%M%S"`
iter=4

for (( i=1; i<=$iter; i++ ))
do
   echo -e "\n\nlaunching $i / $iter\n"
   ./ising -o ${mydir}_${i}
done

