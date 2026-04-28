#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=3:00:00
#SBATCH --job-name=protgMTD
#SBATCH --output=../Results/Rates/rates_metad.out

colvar="metad_short.colvar"
log_file="p.log"
paces=(1e4 2e4 5e4 1e5 5e5 1e6)

additional_options="--temp 312 --numboots 1000 --seed 54321 -bq --nooffset"

head_dir="../../Data/protG/Q_frac_native_contacts_metad/qruns_pace"
command=""
for pace in "${paces[@]}"
do
    command+=" -i ${head_dir}${pace}/run_*/${colvar} --barrier $pace --logfiles ${head_dir}${pace}/run_*/${log_file}"
done
echo $head_dir
python ../Software/rates_eatr_opes.py $command --vcol 2 $additional_options


colvar="metad.colvar"
paces=(1e4 2e4 5e4 1e5 5e5 1e6)

head_dir="../../Data/protG/E_end_end_distance_metad/eruns_pace"
command=""
for pace in "${paces[@]}"
do
    command+=" -i ${head_dir}${pace}/run_*/${colvar} --barrier $pace --logfiles ${head_dir}${pace}/run_*/${log_file}"
done
echo $head_dir
python ../Software/rates_eatr_opes.py $command --vcol 2 $additional_options
