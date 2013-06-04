#!/bin/bash

intdir=L12H0_XXX
outdir=xxx
status=("q" "qk1r" "qk1i" "qk2r" "qk2i")

i=0
while [ $i -lt ${#status[@]} ]
do
    c=0
    for file in ./$(indir}*/${status[i]}_*
    do
	cc=`printf "%04d" $c`
	cp $file ${outdir}/${status[i]}_pair_${cc}.txt
	#echo "cp $file ${outdir}/${status[i]}_pair_${cc}.txt"
	let c=$(($c + 1))
    done

    let i=$(($i + 1))
done


