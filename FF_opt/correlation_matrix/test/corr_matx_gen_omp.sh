#!/bin/bash

complex=1abeA2
echo $complex

#track=${complex}.high_bk
track=1abeA2.high_bk
corr_mat_ofn=${track}.mat_omp
total_row=$(wc -l ${track}|awk '{print $1}')
pdb=$(echo ${complex:0:5}).pdb
ff=${complex}-0.4.ff
sdf=${complex}.sdf

working_dir=./
ff_dir=/work/jaydy/dat/astex/${complex}/

parameter_dir=~/dat/
bin=/home/jaydy/bin/corr_matx_gen_omp

mkdir -p ${working_dir}

cp ${ff_dir}{${pdb},${ff},${sdf}} ${working_dir}
cp ${parameter_dir}{08ff_opt,gpudocksm.ff,nor_a,nor_b,baye_nor_a,baye_nor_b} $working_dir

# running
cd $working_dir

echo "running correlation matrix generation"

temp=0.0044f

rm -rf output_*

${bin} -floor_temp ${temp}f -ceiling_temp 0.036f -nt 1 -t 0.02f -r 0.08f -ns 20000 -nc 1 -p ${pdb} -l ${sdf} -s ${ff} -id ${complex} -tr ${track} -total_row ${total_row} -o ${corr_mat_ofn}


##mark the finish time
date

##and we are out'a here
exit 0
