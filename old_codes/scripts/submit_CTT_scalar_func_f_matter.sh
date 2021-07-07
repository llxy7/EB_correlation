#!/usr/bin/bash

#span=$4

# f_14.96_L_-10.36_tau_1_N_100_kmax_6_nat/ClTT_rei_tensor/func_f/
mkdir -p /tomerv3/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_scalar_matter/func_f

#for ll in {2..10}
#for ll in {20..100..10}
#for ll in {800..1000..100}
#for ll in 2 3 5 8 10 20 30 50 70 100
#for ll in 200
#for ll in 300 500
#for ll in 800 2000
#for ll in 2 5 8 10 30 50 80 100
for ll in 300 500
#for ll in 800 1000 2000
do
	for((kk=1;kk<=$4;kk++))
	do
		qsub -q tomerv -N ll_${ll}_kk_${kk} /tomerv3/Sida/b-mode-log-griding/code/old_codes/scripts/script_CTT_scalar_func_f_matter.PBS -v "L=$2, tau=$3, n=$4, k=$5, ll=${ll}, kk=${kk}" -e /tomerv3/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_scalar_matter/func_f/ll_"${ll}"_kk_${kk}.e -o /tomerv3/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_scalar_matter/func_f/ll_"${ll}"_kk_${kk}.o
	done
done

