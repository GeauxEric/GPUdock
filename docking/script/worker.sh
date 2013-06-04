#!/bin/sh

SET_ID1=$1
SET_TMP=$2

SET_ID2=`echo $SET_TMP | sed -e s/\-/\ /g | awk '{print $1}'`
SET_ID3=`echo $SET_TMP | sed -e s/\-/\ /g | awk '{print $2}'`

export PATH=$PATH:/home/michal/local/bin

export DAT_PATH=/work/michal/antybole

export BIN_PATH=/home/michal/local

TMP_DIR=/var/scratch/$USER/tmp-$SET_ID1-$SET_ID2-$RANDOM

mkdir -p $TMP_DIR

cd $TMP_DIR/

tar -xzf $DAT_PATH/subsets-inp/subset$SET_ID1/part$SET_ID2.tar.gz -C $TMP_DIR/

gunzip -c part$SET_ID2/file$SET_ID3.sdf.gz > file$SET_ID3.sdf

# generate rich sdf

export EFSBABEL=$BIN_PATH/bin/babel

export EFSOBPROP=$BIN_PATH/bin/obprop

export EFSMAYACP=$BIN_PATH/src/mayachemtools/bin/CalculatePhysicochemicalProperties.pl

export EFSMAYAMF=$BIN_PATH/src/mayachemtools/bin/MACCSKeysFingerprints.pl

$BIN_PATH/bin/perl $BIN_PATH/bin/efindsite_sdf -s file$SET_ID3.sdf -o file$SET_ID3-tmp1.sdf -i MOLID -f -c

# construct ensemble

export ESDBABEL=$BIN_PATH/bin/babel

export ESDBALLOON=$BIN_PATH/src/balloon

export ESDSCLUSTER=$BIN_PATH/src/cluto-2.1.2/Linux-x86_64/scluster

$BIN_PATH/bin/perl $BIN_PATH/bin/esimdock_ens -s file$SET_ID3-tmp1.sdf -o file$SET_ID3-tmp2.sdf -i MOLID -n 50

# get output

rm file$SET_ID3-tmp1.sdf file$SET_ID3.sdf

mv file$SET_ID3-tmp2.sdf file$SET_ID3.sdf

gzip -9 file$SET_ID3.sdf

mkdir -p $DAT_PATH/subsets-out/subset$SET_ID1/part$SET_ID2

mv file$SET_ID3.sdf.gz $DAT_PATH/subsets-out/subset$SET_ID1/part$SET_ID2/

cd $TMP_DIR/../

rm -rf $TMP_DIR

exit 0
