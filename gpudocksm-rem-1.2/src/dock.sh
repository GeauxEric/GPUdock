#!/bin/sh
rm -rf output_*

# floor_temp
# lowest temp
# ceiling_temp
# highest temp
# num_temp
# total number of temperatures
# t
# translational scale
# r
# rotation scale


./dock -floor_temp 0.01f -ceiling_temp 0.3f -num_temp 20 -t 0.01f -r 0.08f -s 30000 -stp_per_exchg 10 
