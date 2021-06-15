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

echo $1 $2 $3 $4 $5

mkdir -p f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClBB_rei_accurate

for((uu=1;uu<=Ntot;uu++))
do
        mkdir -p f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClBB_rei_accurate/uu_${uu}_GF 
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
               	#tmp=$((uu-vv>0 ? uu-vv : vv-uu));
               	#for((ww=$((tmp>1 ? tmp : 1));ww<=$((uu+vv<Ntot ? uu+vv : Ntot));ww++))
               	#do
		qsub -q tomerv -N uu_${uu}_vv_${vv} /tomerv/tomerv-shared/Sida/b-mode-log-griding/script_CBB_rei_accurate_GF.PBS -v "f=$1, L=$2, tau=$3, n=$4, k=$5, uu=${uu}, vv=${vv}" -e /tomerv/tomerv-shared/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClBB_rei_accurate/uu_${uu}_GF/uu_${uu}_vv_${vv}.e -o /tomerv/tomerv-shared/Sida/b-mode-log-griding/f_"$1"_L_"$2"_tau_"$3"_N_"$4"_kmax_"$5"_nat_tau/ClBB_rei_accurate/uu_${uu}_GF/uu_${uu}_vv_${vv}.o
                #done
	done
	((Nsubmitted+=Ntosubmit))
	sleep 1m
done

#qstat -u SENSEI1 > log_queue_status.log
#Njobs=$( wc -l < log_queue_status.log )

#while((Njobs>0))
#do
#        sleep 1m
#        qstat -u SENSEI1 > log_queue_status.log
#        Njobs=$( wc -l < log_queue_status.log )
#done

#./collect.sh $1 $2 $3 $4 $5 > log_collect.log

#while [ -f ERROR.err ]
#do
#	./resubmit.sh $1 $2 $3 $4 $5 > log_resubmit.log
#	qstat -u SENSEI1 > log_queue_status.log
#	Njobs=$( wc -l < log_queue_status.log )
#
#	while((Njobs>0))
#	do
#        	sleep 1m
#        	qstat -u SENSEI1 > log_queue_status.log
#        	Njobs=$( wc -l < log_queue_status.log )
#	done
#	./collect.sh $1 $2 $3 $4 $5 > log_collect.log
#done

#math -script sum.m  $1 $2 $3 $4 $4 > log_sum.log

