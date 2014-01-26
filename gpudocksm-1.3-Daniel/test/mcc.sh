
pdb=1atlB.pdb
sdf=1atlB1.sdf
matrix_file=1atlB1-A-matrix
decoy_file=1atlB1-A-Decoys.txt

exe=/usr/local/packages/mpich/3.0.2/Intel-13.0.0/bin/mpiexec.hydra

echo "$complex mcc matrix generating ... begin"

export GPUDOCKSMDAT_FF=$HOME/work/dat/gpudocksm.ff
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/michal/apps/gsl/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/packages/openmpi/1.6.2/gcc-4.4.6/lib
export NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`

# $exe -f $PBS_NODEFILE -n $NPROCS ./mpi_mcc -p $pdb -l $sdf -o ${matrix_file}1 -d ${decoy_file}1
$exe -f $PBS_NODEFILE -n 2 ./mpi_mcc -p $pdb -l $sdf -o ${matrix_file} -d ${decoy_file}
