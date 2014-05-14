#!/bin/bash

##PBS -q gpu
##PBS -A hpc_docking
#PBS -q lasigma

#PBS -j oe

## user specified
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#send an email after job aborts 
#PBS -M yding8@lsu.edu
#PBS -m a


date

# complex=$var1
# complex=1dwbB1
# complex=1abeA2
# complex=1abfA1
complex=1a07C1
# complex=$var1
echo $complex

pdb=$(echo ${complex:0:5}).pdb
ff=${complex}-0.8.ff
sdf=${complex}.sdf
opt_fn=08ff_opt
nor_a_fn=08_nor_a
nor_b_fn=08_nor_b

# working_dir=~/work/working/linear_z_score/${complex}/
# working_dir=/work/jaydy/working/08ff_z_score/${complex}/
working_dir=/home/jaydy/src/GPUdock/gpudocksm-rem-1.2/test/
ff_dir=/work/jaydy/dat/astex/${complex}/

parameter_dir=/work/jaydy/dat/
bin=/home/jaydy/src/GPUdock/gpudocksm-rem-1.2/test/linear_dock

mkdir -p ${working_dir}

cp ${ff_dir}{${pdb},${ff},${sdf}} ${working_dir}
cp ${parameter_dir}{${opt_fn},gpudocksm.ff,${nor_a_fn},${nor_b_fn},baye_nor_a,baye_nor_b} $working_dir

# running
cd $working_dir

echo "generate ligand conf on the fly"

# temp=0.044f
temp=0.044f

rm -rf output_*

# ${bin} -floor_temp ${temp}f -ceiling_temp 0.036f -nt 1 -t 0.02f -r 0.08f -ns 3000 -nc 1 -p ${pdb} -l ${sdf} -s ${ff} -id ${complex} -opt_fn ${opt_fn} > report
${bin} -floor_temp ${temp}f -ceiling_temp 0.036f -nt 1 -t 0.02f -r 0.08f -ns 3000 -nc 1 -p ${pdb} -l ${sdf} -s ${ff} -id ${complex} -opt_fn ${opt_fn} -tr ligand_pos.in


##mark the finish time
date

##and we are out'a here
exit 0
