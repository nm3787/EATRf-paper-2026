#!/bin/bash

# Create the STATE files for the biased PMFs figure.

for barr in 5 7 9 11 13
do
	for run in {1..100}
	do
		dir_prefix="../Data/protG/Q_frac_native_contacts_opes/qruns_barr${barr}/run_${run}"
		python State_from_Kernels.py --kernels $dir_prefix/KERNELS --outfile $dir_prefix/STATE
		dir_prefix="../Data/protG/E_end_end_distance_opes/eruns_barr${barr}/run_${run}"
                python State_from_Kernels.py --kernels $dir_prefix/KERNELS --outfile $dir_prefix/STATE
	done
done
