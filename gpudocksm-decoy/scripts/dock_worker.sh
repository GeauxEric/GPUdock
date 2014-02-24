#!/bin/sh

SET_ID1=$1
SET_TMP=$2
export BIN_PATH=$3


SET_ID2=`echo $SET_TMP | sed -e s/\-/\ /g | awk '{print $1}'`
SET_ID4=`echo $SET_ID2 | cut -c1-5`

##export PATH=$PATH:../bin
export DAT_PATH=/work/$USER

TMP_DIR=/var/scratch/$USER/tmp-$SET_ID2-$RANDOM

mkdir -p $TMP_DIR

cd $TMP_DIR/

cp $BIN_PATH/dat/gpudocksm.ff $TMP_DIR/

export GPUDOCKSMDAT_FF=gpudocksm.ff
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/michal/apps/gsl/lib/

cp $DAT_PATH/astex/$SET_ID2/$SET_ID4.pdb $TMP_DIR/
cp $DAT_PATH/astex/$SET_ID2/$SET_ID2.sdf $TMP_DIR/
cp $DAT_PATH/astex/$SET_ID2/$SET_ID2-0.4.ff $TMP_DIR/
##tar -xzf $DAT_PATH/proteins/$SET_ID1.tar.gz -C $TMP_DIR/

##gunzip -c $SET_ID1.gz > file$SET_ID3

##$BIN_PATH/bin/gpudocksm -p $pdb -l $sdf -s $ff -i MOLID -ns 10 -nc 10
$BIN_PATH/bin/gpudocksm -p $SET_ID4.pdb -l $SET_ID2.sdf -s $SET_ID2-0.4.ff -i MOLID -ns 130 -nc 130

mv decoyeneH.txt $SET_ID2-eneH.txt
mv decoyeneL1.txt $SET_ID2-eneL1.txt
mv decoyeneL2.txt $SET_ID2-eneL2.txt
mv decoysH.txt $SET_ID2-decoysH.txt
mv decoysL1.txt $SET_ID2-decoysL1.txt
mv decoysL2.txt $SET_ID2-decoysL2.txt
mv MCcoefsH.txt $SET_ID2-mccH.txt
mv MCcoefsL1.txt $SET_ID2-mccL1.txt
mv MCcoefsL2.txt $SET_ID2-mccL2.txt

gzip -9 $SET_ID2-eneH.txt
gzip -9 $SET_ID2-eneL1.txt
gzip -9 $SET_ID2-eneL2.txt
gzip -9 $SET_ID2-decoysH.txt
gzip -9 $SET_ID2-decoysL1.txt
gzip -9 $SET_ID2-decoysL2.txt
gzip -9 $SET_ID2-mccH.txt
gzip -9 $SET_ID2-mccL1.txt
gzip -9 $SET_ID2-mccL2.txt

mkdir -p $DAT_PATH/decoys/$SET_ID2

mv $SET_ID2-eneH.txt.gz $DAT_PATH/decoys/$SET_ID2/
mv $SET_ID2-eneL1.txt.gz $DAT_PATH/decoys/$SET_ID2/
mv $SET_ID2-eneL2.txt.gz $DAT_PATH/decoys/$SET_ID2/
mv $SET_ID2-decoysH.txt.gz $DAT_PATH/decoys/$SET_ID2/
mv $SET_ID2-decoysL1.txt.gz $DAT_PATH/decoys/$SET_ID2/
mv $SET_ID2-decoysL2.txt.gz $DAT_PATH/decoys/$SET_ID2/
mv $SET_ID2-mccH.txt.gz $DAT_PATH/decoys/$SET_ID2/
mv $SET_ID2-mccL1.txt.gz $DAT_PATH/decoys/$SET_ID2/
mv $SET_ID2-mccL2.txt.gz $DAT_PATH/decoys/$SET_ID2/


cd $TMP_DIR/../

rm -rf $TMP_DIR

exit 0
