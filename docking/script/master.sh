#!/bin/sh

export PATH=$PATH:/home/michal/local/bin

export DAT_PATH=/work/michal/antybole

export BIN_PATH=/home/michal/local

export SET1=000
export SET2=01

cd $DAT_PATH/

cat $PBS_NODEFILE | grep qb | sort | uniq | awk '{print "8/"$1}' > machines-$SET1-$SET2.lst

# rich sdf + ensemble

cat subsets-inp/subset$SET1-$SET2.lst | parallel --sshloginfile machines-$SET1-$SET2.lst sh $DAT_PATH/scripts/worker.sh $SET1 {} > /dev/null 2>&1

rm machines-$SET1-$SET2.lst

exit 0
