#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=3:00:00
#SBATCH --job-name=protGopes
#SBATCH --output=../Results/Rates/rates_protg_opes.out

colvar="opes_short.colvar"
log_file="p.log"
barrs=(5 7 9 11 13)

additional_options="--temp 312 --numboots 1000 --seed 12345 --opesf --cdf -bq"

head_dir="../Data/protG/Q_frac_native_contacts_opes/qruns_barr"
command=""
for barr in "${barrs[@]}"
do
    command+=" -i ${head_dir}${barr}/run_*/${colvar} --barrier $barr --logfiles ${head_dir}${barr}/run_*/${log_file}"
done
echo $head_dir
python ../Software/rates_eatr_opes.py $command --vcol 3 $additional_options

head_dir="../Data/protG/E_end_end_distance_opes/eruns_barr"
command=""
for barr in "${barrs[@]}"
do
    command+=" -i ${head_dir}${barr}/run_*/${colvar} --barrier $barr --logfiles ${head_dir}${barr}/run_*/${log_file}"
done
echo $head_dir
python ../Software/rates_eatr_opes.py $command --vcol 4 $additional_options

head_dir="../Data/protG/G_radius_of_gyration_opes/gruns_barr"
command=""
for barr in "${barrs[@]}"
do
    command+=" -i ${head_dir}${barr}/run_*/${colvar} --barrier $barr --logfiles ${head_dir}${barr}/run_*/${log_file}"
done
echo $head_dir
python ../Software/rates_eatr_opes.py $command --vcol 4 $additional_options

head_dir="../Data/protG/2D_GE_opes/geruns_barr"
command=""
for barr in "${barrs[@]}"
do
    command+=" -i ${head_dir}${barr}/run_*/${colvar} --barrier $barr --logfiles ${head_dir}${barr}/run_*/${log_file}"
done
echo $head_dir
python ../Software/rates_eatr_opes.py $command --vcol 5 $additional_options
