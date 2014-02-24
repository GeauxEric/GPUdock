export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/michal/apps/gsl/lib/
export GPUDOCKSMDAT_FF=$HOME/work/dat/gpudocksm.ff


# for testing 
# cd ~/work/test/
# ns is monte carlo steps 
# nc is replica exchange cycles
# t is temperature
# nr is number of replicas
# bz is the boltzman constant

$HOME/work/bin/decoy_gen -p 1a07C.pdb -l 1a07C1.sdf -s 1a07C1-0.4.ff -i MOLID -ns 10 -nc 10 -nr 6 -t 1.1 -bz 0.7 -dp 1000


##############################
# for real simulation 
# echo "decoys generating starting..."
# 
# prt=$1
# complex=$2 
# pdb=$prt.pdb
# sdf=$complex.sdf
# ff=$3
# mc_stp=$4
# rep_cycles=$5
# tmp_dir=$HOME/work/working/$complex
# 
# # jump to tmp_dir, preparing to execute the program
# cd $tmp_dir
# echo "cd $tmp_dir"
# 
# 
# $HOME/work/bin/decoy_gen -p $pdb -l $sdf -s $ff -i MOLID -ns $mc_stp -nc $rep_cycles -nr 6 -t 1.1 -bz 0.6
