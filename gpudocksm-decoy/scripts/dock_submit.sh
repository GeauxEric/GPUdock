
qsub -q workq -A lasigma -l walltime=48:00:00 -o gpudock-01.log -e gpudock-01.log -j oe -l nodes=1:ppn=16 dock_master1.sh
