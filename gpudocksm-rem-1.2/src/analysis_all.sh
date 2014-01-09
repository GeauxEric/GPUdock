################################################################################
# use ./analysis to load the hdf5 record files and redirect the output to a csv format
num_temp=$(cat report | grep "total temperatures"|awk '{print $3}')
num_lig_conf=$(cat report |grep "total ligand conformations"|awk '{print $4}')
num_prt_conf=$(cat report |grep "total prt conformations"|awk '{print $4}')
temps=$(cat report|grep "temp #"|awk '{print $4}')

prt_conf=0
lig_conf=0
i=0

for temp in $temps
do
    let "rep_num = $num_temp * $num_lig_conf * $prt_conf + $num_lig_conf * $i + $lig_conf"
    printf "analyze replica # %d with temperature %f\n" "$rep_num" "$temp"
    sh analysis.sh $rep_num $temp

    i=$((i+1))
done

# for ((i=0; i < $num_temp; i++))
# do
#     # to calculate flatten addr of replica
#     let "num_rep = $num_temp * $num_lig_conf * $prt_conf + $num_lig_conf * $i + $lig_conf"
#     # call analysis.sh to analyze that replica
#     sh analysis.sh $num_rep
# done
