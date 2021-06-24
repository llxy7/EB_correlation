#!/usr/bin/bash

#span=$4

# f_14.96_L_-10.36_tau_1_N_100_kmax_6_nat/ClTT_rei_tensor/func_f/
mkdir -p /tomerv3/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/GF_new

for((kk=1;kk<=$4;kk++))
do
	qsub -q tomerv -N kk_${kk} /tomerv3/Sida/b-mode-log-griding/code/old_codes/scripts/script_GF_tensor.PBS -v "f=$1, L=$2, tau=$3, n=$4, k=$5, kk=${kk}" -e /tomerv3/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/GF_new/kk_${kk}.e -o /tomerv3/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/GF_new/kk_${kk}.o
done

