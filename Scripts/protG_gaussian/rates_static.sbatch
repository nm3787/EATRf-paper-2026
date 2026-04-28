#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=03:00:00
#SBATCH --job-name=protGgauss
#SBATCH --output=../../Results/Rates/rates_static.out


echo "Q"
head_dir="../../Data/protG/Q_frac_native_contacts_gaussian/qruns_"
python rates_eatr_static.py -i ${head_dir}4kJ/*/bias.colvar --barrier 4 --logfiles ${head_dir}4kJ/*/p.log -i ${head_dir}5kJ/*/bias.colvar --barrier 5 --logfiles ${head_dir}5kJ/*/p.log -i ${head_dir}6kJ/*/bias.colvar --barrier 6 --logfiles ${head_dir}6kJ/*/p.log -i ${head_dir}7kJ/*/bias.colvar --barrier 7 --logfiles ${head_dir}7kJ/*/p.log -i ${head_dir}8kJ/*/bias.colvar --barrier 8 --logfiles ${head_dir}8kJ/*/p.log -i ${head_dir}10kJ/*/bias.colvar --barrier 10 --logfiles ${head_dir}10kJ/*/p.log --temp 312 --tcol 0 --vcol 4 -bq --numboots 1000 --seed 9586 --opesf --nooffset --cdf --cv Q

echo "E"
head_dir="../../Data/protG/E_end_end_distance_gaussian/eruns_"
python rates_eatr_static.py -i ${head_dir}4kJ/*/bias.colvar --barrier 4 --logfiles ${head_dir}4kJ/*/p.log -i ${head_dir}5kJ/*/bias.colvar --barrier 5 --logfiles ${head_dir}5kJ/*/p.log -i ${head_dir}6kJ/*/bias.colvar --barrier 6 --logfiles ${head_dir}6kJ/*/p.log -i ${head_dir}7kJ/*/bias.colvar --barrier 7 --logfiles ${head_dir}7kJ/*/p.log -i ${head_dir}8kJ/*/bias.colvar --barrier 8 --logfiles ${head_dir}8kJ/*/p.log -i ${head_dir}10kJ/*/bias.colvar --barrier 10 --logfiles ${head_dir}10kJ/*/p.log --temp 312 --tcol 0 --vcol 4 -bq --numboots 1000 --seed 5439 --opesf --nooffset --cdf --cv E
