#!/bin/sh

# # floor_temp
# # lowest temp
# # ceiling_temp
# # highest temp
# # num_temp
# # total number of temperatures
# # t
# # translational scale
# # r
# # rotation scale


printf "%s\t\t%s\n" "stp_per_exchg" "exchg_AR"

for i in {10..200..10}
do
    printf "%d\t\t" "$i"

    rm -rf output_*  # remove older data
    ./dock -floor_temp 0.1f -ceiling_temp 0.3f -num_temp 2 -t 0.01f -r 0.08f -s 60000 -stp_per_exchg $i \
	|grep "AR of temp exchange"|awk '{print $8}'
done 
