export GPUDOCKSMDAT_FF=$HOME/work/dat/gpudocksm.ff
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/michal/apps/gsl/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/packages/openmpi/1.6.2/gcc-4.4.6/lib/

##########################
# for testing
# cd $HOME/work/test/
# export NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`
# mpirun -machinefile $PBS_NODEFILE -np $NPROCS $HOME/work/bin/mpi_mcc -p 1a07C.pdb -l 1a07C1.sdf -o 1a07C1-A-matrix.txt -d 1a07C1-A-Decoys.txt
# mpirun -np 16 $HOME/work/bin/mpi_mcc -p 1a07C.pdb -l 1a07C1.sdf -o 1a07C1-A-matrix.txt0 -d 1a07C1-A-Decoys.txt0
# mpirun -np 16 $HOME/work/bin/mpi_mcc -p 1a07C.pdb -l 1a07C1.sdf -o 1a07C1-A-matrix.txt1 -d 1a07C1-A-Decoys.txt1
# mpirun -np 16 $HOME/work/bin/mpi_mcc -p 1a07C.pdb -l 1a07C1.sdf -o 1a07C1-A-matrix.txt2 -d 1a07C1-A-Decoys.txt2
# mpirun -np 16 $HOME/work/bin/mpi_mcc -p 1a07C.pdb -l 1a07C1.sdf -o 1a07C1-A-matrix.txt3 -d 1a07C1-A-Decoys.txt3

mpi_run=/usr/local/packages/mpich/3.0.2/Intel-13.0.0/bin/mpirun
cd $HOME/work/test/

export NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`
$mpi_run -f $PBS_NODEFILE -n $NPROCS $HOME/work/bin/mpi_mcc_test -p 1a07C.pdb -l 1a07C1.sdf -o 1a07C1-A-matrix.txt -d 1a07C1-A-Decoys.txt

#########################
# for real simulation

# echo "mcc matrix generating ..."
# complex=$1
# prt=$2
# matrix_file=$3
# decoy_file=$4 
# pdb=$prt.pdb
# sdf=$complex.sdf
# tmp_dir=$HOME/work/working/$complex
# # tmp_dir=$HOME/work/test/
# cp $HOME/work/bin/mpi_mcc $tmp_dir
# 
# export NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`
# 
# 
# # jump to tmp_dir, preparing to execute the program
# cd $tmp_dir
# echo "cd $tmp_dir"
# mpirun -np $NPROCS ./mpi_mcc -p $pdb -l $sdf -o $matrix_file -d $decoy_file

# ##################
# for mpi testing
# cd ~/work/bin/
# 
# export NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`
# mpirun -np $NPROCS ./mpi_mcc
