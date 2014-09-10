#!/bin/bash

#complex=$var1
#echo $complex

complex=1a07C1

bin=../src/dock
input_dir=.
parameters_dir=${input_dir}/parameters
complex_dir=${input_dir}/${complex}

pdb_file=${complex_dir}/$(echo ${complex:0:5}).pdb
sdf_file=${complex_dir}/${complex}.sdf
ff_file=${complex_dir}/${complex}-0.8.ff

opt_file=${parameters_dir}/08ff_opt
nor_a_file=${parameters_dir}/08_nor_a
nor_b_file=${parameters_dir}/08_nor_b
para_file=${parameters_dir}/paras

near_native_trace_file=/work/jaydy/GitHub/GPUdock/SimilarityMatirxGen/data/test_matrix.txt
near_native_simi_file=/work/jaydy/GitHub/GPUdock/SimilarityMatirxGen/data/test_simi_matrix.txt


cmd="\
${bin} \
-id ${complex} \
-p ${pdb_file} \
-l ${sdf_file} \
-s ${ff_file} \
\
-opt ${opt_file} \
-na ${nor_a_file} \
-nb ${nor_b_file} \
-para ${para_file} \
\
-ns 1000 \
-nc 10 \
-floor_temp 0.044f \
-ceiling_temp 0.036f \
-nt 1 \
-t 0.02f \
-r 0.08f \
\
-ntr ${near_native_trace_file} \
-simi ${near_native_simi_file}
"


rm -rf out*
echo ${cmd}
${cmd}




