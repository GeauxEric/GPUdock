#!bin/bash
min_temp=0.0046080000
step=1.1
temp=$min_temp


printf "temperature\t\tAR_MC\n"
for ((i=1; i < 20; i += 1))
do
    printf "%.10f\t\t" "$temp"
    ./bayes_dock -floor_temp ${temp}f -ceiling_temp 0.036f -nt 1 -t 0.01f -r 0.04f -ns 300000 -nc 1 -p 1a07C.pdb -l 1a07C1.sdf -s 1a07C1.ff -id 1a07C1 |grep "AR of MC" |awk '{print $7}'

    temp=$(echo "$temp * $step"|bc)
done
