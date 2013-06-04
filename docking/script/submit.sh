SUBSET=000

qsub -q workq -A loni_csbg02 -M mbrylinski@lsu.edu -l walltime=72:00:00 -o s02-$SUBSET-01.log -e s02-$SUBSET-01.log -j oe -l nodes=5:ppn=8 s02-$SUBSET-01
qsub -q workq -A loni_csbg02 -M mbrylinski@lsu.edu -l walltime=72:00:00 -o s02-$SUBSET-02.log -e s02-$SUBSET-02.log -j oe -l nodes=5:ppn=8 s02-$SUBSET-02
qsub -q workq -A loni_csbg02 -M mbrylinski@lsu.edu -l walltime=72:00:00 -o s02-$SUBSET-03.log -e s02-$SUBSET-03.log -j oe -l nodes=5:ppn=8 s02-$SUBSET-03
qsub -q workq -A loni_csbg02 -M mbrylinski@lsu.edu -l walltime=72:00:00 -o s02-$SUBSET-04.log -e s02-$SUBSET-04.log -j oe -l nodes=5:ppn=8 s02-$SUBSET-04
qsub -q workq -A loni_csbg02 -M mbrylinski@lsu.edu -l walltime=72:00:00 -o s02-$SUBSET-05.log -e s02-$SUBSET-05.log -j oe -l nodes=5:ppn=8 s02-$SUBSET-05
qsub -q workq -A loni_csbg02 -M mbrylinski@lsu.edu -l walltime=72:00:00 -o s02-$SUBSET-06.log -e s02-$SUBSET-06.log -j oe -l nodes=5:ppn=8 s02-$SUBSET-06
qsub -q workq -A loni_csbg02 -M mbrylinski@lsu.edu -l walltime=72:00:00 -o s02-$SUBSET-07.log -e s02-$SUBSET-07.log -j oe -l nodes=5:ppn=8 s02-$SUBSET-07
qsub -q workq -A loni_csbg02 -M mbrylinski@lsu.edu -l walltime=72:00:00 -o s02-$SUBSET-08.log -e s02-$SUBSET-08.log -j oe -l nodes=5:ppn=8 s02-$SUBSET-08
qsub -q workq -A loni_csbg02 -M mbrylinski@lsu.edu -l walltime=72:00:00 -o s02-$SUBSET-09.log -e s02-$SUBSET-09.log -j oe -l nodes=5:ppn=8 s02-$SUBSET-09
qsub -q workq -A loni_csbg02 -M mbrylinski@lsu.edu -l walltime=72:00:00 -o s02-$SUBSET-10.log -e s02-$SUBSET-10.log -j oe -l nodes=5:ppn=8 s02-$SUBSET-10
