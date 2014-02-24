#!/bin/bash
## To run on SMII
## There are 204 proteins to run, including ones that include a metal atom
#PBS -q workq
#PBS -A lasigma
#
#PBS -l nodes=1:ppn=16
#
#PBS -l walltime=48:00:00
#
#PBS -o /scratch/dcase1/docking/output
#PBS -j oe
#
#PBS -N dcaseDocking5
date

##export HOME_DIR=/home/dcase1/docking
export WORK_DIR=/home/dcase1/biology/trunk/gpudocksm-1.3.1-Daniel/scripts

##cp $HOME_DIR/hydro $WORK_DIR
cd $WORK_DIR
echo $PWD
bash dock_master5.sh

date

exit 0
