#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --job-name=cavopes
#SBATCH --output=../Results/Rates/rates_cavity_opes.out


echo 0.87z
head_dir="../Data/cavity/Z_87_cavity_opes"
python ../Software/rates_eatr_opes.py -i $head_dir/barrier3.5/*.opes.out --barrier 3.5 --logfiles $head_dir/barrier3.5/*.log -i $head_dir/barrier4/*.opes.out --barrier 4 --logfiles $head_dir/barrier4/*.log -i $head_dir/barrier4.5/*.opes.out --barrier 4.5 --logfiles $head_dir/barrier4.5/*.log -i $head_dir/barrier5/*.opes.out --barrier 5 --logfiles $head_dir/barrier5/*.log -i $head_dir/barrier5.5/*.opes.out --barrier 5.5 --logfiles $head_dir/barrier5.5/*.log -i $head_dir/barrier6/*.opes.out --barrier 6 --logfiles $head_dir/barrier6/*.log --beta 1.6773984449958872 --tcol 0 --vcol 1 --timeunit 1e-15 -bq --numboots 1000 --opesf --seed 12345

echo 0.9z
head_dir="../Data/cavity/Z_90_cavity_opes"
python ../Software/rates_eatr_opes.py -i $head_dir/barrier3.5/*.opes.out --barrier 3.5 --logfiles $head_dir/barrier3.5/*.log -i $head_dir/barrier4/*.opes.out --barrier 4 --logfiles $head_dir/barrier4/*.log -i $head_dir/barrier4.5/*.opes.out --barrier 4.5 --logfiles $head_dir/barrier4.5/*.log -i $head_dir/barrier5/*.opes.out --barrier 5 --logfiles $head_dir/barrier5/*.log -i $head_dir/barrier5.5/*.opes.out --barrier 5.5 --logfiles $head_dir/barrier5.5/*.log -i $head_dir/barrier6/*.opes.out --barrier 6 --logfiles $head_dir/barrier6/*.log --beta 1.6773984449958872 --tcol 0 --vcol 1 --timeunit 1e-15 -bq --numboots 1000 --opesf --seed 23451

echo 0.95z
head_dir="../Data/cavity/Z_95_cavity_opes"
python ../Software/rates_eatr_opes.py -i $head_dir/barrier3.5/*.opes.out --barrier 3.5 --logfiles $head_dir/barrier3.5/*.log -i $head_dir/barrier4/*.opes.out --barrier 4 --logfiles $head_dir/barrier4/*.log -i $head_dir/barrier4.5/*.opes.out --barrier 4.5 --logfiles $head_dir/barrier4.5/*.log -i $head_dir/barrier5/*.opes.out --barrier 5 --logfiles $head_dir/barrier5/*.log -i $head_dir/barrier5.5/*.opes.out --barrier 5.5 --logfiles $head_dir/barrier5.5/*.log -i $head_dir/barrier6/*.opes.out --barrier 6 --logfiles $head_dir/barrier6/*.log --beta 1.6773984449958872 --tcol 0 --vcol 1 --timeunit 1e-15 -bq --numboots 1000 --opesf --seed 34512

echo z
head_dir="../Data/cavity/Z_100_cavity_opes"
python ../Software/rates_eatr_opes.py -i $head_dir/barrier3/*.opes.out --barrier 3 --logfiles $head_dir/barrier3/*.log -i $head_dir/barrier3.5/*.opes.out --barrier 3.5 --logfiles $head_dir/barrier3.5/*.log -i $head_dir/barrier4/*.opes.out --barrier 4 --logfiles $head_dir/barrier4/*.log -i $head_dir/barrier4.5/*.opes.out --barrier 4.5 --logfiles $head_dir/barrier4.5/*.log -i $head_dir/barrier5/*.opes.out --barrier 5 --logfiles $head_dir/barrier5/*.log -i $head_dir/barrier6/*.opes.out --barrier 6 --logfiles $head_dir/barrier6/*.log --beta 1.6773984449958872 --tcol 0 --vcol 1 --timeunit 1e-15 -bq --numboots 1000 --opesf --seed 45123

echo rho
head_dir="../Data/cavity/P_rho_cavity_opes"
python ../Software/rates_eatr_opes.py -i $head_dir/barrier3/*.opes.out --barrier 3 --logfiles $head_dir/barrier3/*.log -i $head_dir/barrier3.5/*.opes.out --barrier 3.5 --logfiles $head_dir/barrier3.5/*.log -i $head_dir/barrier4/*.opes.out --barrier 4 --logfiles $head_dir/barrier4/*.log -i $head_dir/barrier4.5/*.opes.out --barrier 4.5 --logfiles $head_dir/barrier4.5/*.log -i $head_dir/barrier5/*.opes.out --barrier 5 --logfiles $head_dir/barrier5/*.log -i $head_dir/barrier6/*.opes.out --barrier 6 --logfiles $head_dir/barrier6/*.log --beta 1.6773984449958872 --tcol 0 --vcol 1 --timeunit 1e-15 -bq --numboots 1000 --opesf --seed 51234
