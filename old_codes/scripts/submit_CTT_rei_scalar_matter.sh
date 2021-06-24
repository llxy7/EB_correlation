#!/usr/bin/bash

if [ -z $6 ]
then
        echo "give me 6 parameters: log10f, log10Lambda, pushback, Ntot, coeffk, JobNumInQueue"
        exit 1
fi

Ntot=$4
#span=$4

MaxInQueue=$(($6+5))

MaxJobNum=$((Ntot*Ntot))

Nsubmitted=0

mkdir -p /tomerv3/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_scalar_matter

for((uu=1;uu<=Ntot;uu++))
do
        mkdir -p /tomerv3/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_scalar_matter/uu_${uu} 
done

while((Nsubmitted<MaxJobNum))
do
	qstat -u sidalu > log_queue_status.log
	Njobs=$( wc -l < log_queue_status.log )
	Ntosubmit=$((MaxInQueue-Njobs))
	cnt=0
	while((cnt<Ntosubmit && ((Nsubmitted+cnt))<MaxJobNum))
	do
		uu=$(((Nsubmitted+cnt)/Ntot+1))
		vv=$(((Nsubmitted+cnt)%Ntot+1))
		((cnt++))
		echo "${uu} ${vv}"
		qsub -q tomerv -N uu_${uu}_vv_${vv} /tomerv3/Sida/b-mode-log-griding/code/old_codes/scripts/script_CTT_rei_scalar_matter.PBS -v "f=$1, L=$2, tau=$3, n=$4, k=$5, uu=${uu}, vv=${vv}" -e /tomerv3/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_scalar_matter/uu_${uu}/uu_${uu}_vv_${vv}.e -o /tomerv3/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClTT_rei_scalar_matter/uu_${uu}/uu_${uu}_vv_${vv}.o
	done
	((Nsubmitted+=Ntosubmit))
	sleep 1m
done

