$MinPrecision=100;
Off[NIntegrate::ncvb,NIntegrate::slwcon,NIntegrate::eincr,General::munfl];

log10f=ToExpression[$CommandLine[[-5]]];
log10Lambda=ToExpression[$CommandLine[[-4]]];
pushback=ToExpression[$CommandLine[[-3]]];
Ntot=ToExpression[$CommandLine[[-2]]];
coeffk=ToExpression[$CommandLine[[-1]]];

workdir="/tomerv/tomerv-shared/Sida/b-mode-log-griding/f_"<>$CommandLine[[-5]]<>"_L_"<>$CommandLine[[-4]]<>"_tau_"<>$CommandLine[[-3]]<>"_N_"<>$CommandLine[[-2]]<>"_kmax_"<>$CommandLine[[-1]]<>"_nat_tau/ClBB_rei_accurate/";
dumpdir="/tomerv/tomerv-shared/Sida/b-mode-log-griding/f_"<>$CommandLine[[-5]]<>"_L_"<>$CommandLine[[-4]]<>"_tau_"<>$CommandLine[[-3]]<>"_N_"<>$CommandLine[[-2]]<>"_kmax_"<>$CommandLine[[-1]]<>"_nat_tau/";

GeVtog=SetPrecision[1.783*10^-24,100];
GeVtom=SetPrecision[1.973*10^-16,100];
G=SetPrecision[6.71*10^-39,100];
Mpl=1/Sqrt[8π G];
H0=SetPrecision[67.8/(299792458/1000)/(3.086*10^22)*GeVtom,100];
α=1;
τ0=1/H0;(*RD*)
τ0=2/H0;(*MD*)
τr=4*10^40;

Λ=10^log10Lambda;
f=10^log10f;
τosc=5*10^40;
kmax=coeffk*1/(2 Pi) Log[9/(Factorial[7]/(2^21 Pi^2) α) (Mpl^4 Λ^4)/(Λ^4/2)^2]*4/τosc;
kmax=3*10^(-39);
τmin=τosc;
tmax=(*τ0/3*)5τosc;
τr=4*10^40;
τrei=2/(3*H0);
tmax=5*τrei;

optdep$rei=0.08;

lmode=Table[i,{i,1,9}]~Join~Table[10i,{i,1,9}]~Join~Table[100i,{i,1,5}];
allnum=Table[True,Length[lmode]];

ans=Table[0,Length[lmode]];

For[uu=1,uu<=Ntot,uu++,
For[vv=1,vv<=Ntot,vv++,
tab=Import[workdir<>"uu_"<>ToString[uu]<>"_GF/uu_"<>ToString[uu]<>"_vv_"<>ToString[vv]<>".dat","Table"]//Quiet;
If[tab==$Failed,Print[uu,vv];Continue[]];
tab=Flatten[tab];
If[NumberQ[#]&/@tab != allnum,Print[ToString[uu]<>" "<>ToString[vv]<>" not number"];Continue[]];
ans+=tab;
]
]
Print[ans]
Export[dumpdir<>"ClBB_rei_accurate_GF.dat",MapThread[{#1,#1*(#1+1)*#2/(2 Pi)}&,{lmode,ans}]]
