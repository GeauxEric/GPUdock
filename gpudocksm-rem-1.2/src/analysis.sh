################################################################################
# setting the paths
hd_path=$(ls output_*/*.h5)
echo -e "hdf5 path\t\t\t"$hd_path

dir_path=$(dirname $hd_path)

hd_fn=${hd_path##*/}

base_name=${hd_fn%.h5}

csv_path=$dir_path/$base_name.csv
echo -e "csv output path\t\t\t"$csv_path

################################################################################
# use ./analysis to load the hdf5 record file and redirect the output to a csv format
./analysis -nl 20 -l 1 -p 2 $hd_path > $csv_path

################################################################################
# use pandas to laod the csv file and remove the duplicates
base_path=$dir_path/$base_name

python removeDup.py $base_path
