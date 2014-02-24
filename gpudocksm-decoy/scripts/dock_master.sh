#!/bin/sh

# the complex name list should be like plist$SET2.lst

export PATH=/project/michal/apps/parallel/bin:$PATH

BIN1_PATH=${PWD}
cd ../
export BIN_PATH=${PWD}
cd $BIN1_PATH

export SET1=000
export SET2=03

##cd $DAT_PATH/

cat $PBS_NODEFILE | grep mike | sort | uniq | awk '{print "8/"$1}' > machines-meta16-PART$SET2.lst

cat $PBS_NODEFILE | grep mike | sort | uniq | awk '{print "gethostip", $1}' > tmp1-PART ; sh tmp1-PART | awk '{print "8/"$2}' >> machines-meta16-PART$SET2.lst ; rm tmp1-PART

##cat plist$SET2.lst | parallel --sshloginfile machines-meta16-PART.lst sh $PWD/dock_worker.sh $SET1 {} $BIN_PATH > /dev/null 2>&1

cat plist$SET2.lst | parallel --retries 50 --sshloginfile machines-meta16-PART$SET2.lst sh $PWD/dock_worker.sh $SET1 {} $BIN_PATH 


rm machines-meta16-PART$SET2.lst

exit 0 
