#!/usr/bin/bash

#span=$4

# f_14.96_L_-10.36_tau_1_N_100_kmax_6_nat/ClTT_rei_tensor/func_f/
mkdir -p f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_tensor/func_f_GF

#for ll in 200 300 400 500
#for ll in 600 700
#for ll in 700 2000
#for ll in 500 1000
#for ll in 2 2000
for ll in 2 5 8 10 30 2000
do
	for((kk=1;kk<=$4;kk++))
	do
		if [ ! -f f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_tensor/func_f_GF/ll_${ll}_kk_${kk}_2.dat ]
		then
		#	cp f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_tensor/func_f_GF/ll_${ll}_kk_${kk}.dat f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_tensor/func_f_GF/ll_${ll}_kk_${kk}_2.dat
			qsub -q tomerv -N ll_${ll}_kk_${kk} /tomerv/tomerv-shared/Sida/b-mode-log-griding/script_CTT_tensor_func_f_GF_MMA_2.PBS -v "L=$2, tau=$3, n=$4, k=$5, ll=${ll}, kk=${kk}" -e /tomerv/tomerv-shared/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_tensor/func_f_GF/ll_${ll}_kk_${kk}_2.e -o /tomerv/tomerv-shared/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_tensor/func_f_GF/ll_${ll}_kk_${kk}_2.o
		fi
	done
done
