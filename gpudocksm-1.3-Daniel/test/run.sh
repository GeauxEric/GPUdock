export LD_INCLUDE_PATH=/usr/local/packages/gsl/1.15/Intel-13.0.0/include:$LD_INCLUDE_PATH
export LD_LIBRARY_PATH=/usr/local/packages/gsl/1.15/Intel-13.0.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/packages/mpich/3.0.2/Intel-13.0.0/lib:$LD_LIBRARY_PATH
export GPUDOCKSMDAT_FF=./gpudocksm.ff



# ../bin/gpudocksm -p 1a07C_2.pdb -l 1a07C1.sdf -s 1a07C1.ff -o test1.sdf -i MOLID
./mpi_mcc -p 1atlB.pdb -l 1atlB1.sdf -o 1atlB1-A-matrix -d 1atlB1-A-Decoys.txt
