$MinPrecision=100;
Off[NIntegrate::ncvb,NIntegrate::slwcon,NIntegrate::eincr,General::munfl];

log10f=ToExpression[$CommandLine[[-5]]];
log10Lambda=ToExpression[$CommandLine[[-4]]];
pushback=ToExpression[$CommandLine[[-3]]];
Ntot=ToExpression[$CommandLine[[-2]]];
coeffk=ToExpression[$CommandLine[[-1]]];

workdir="/tomerv3/Sida/b-mode-log-griding/f_"<>$CommandLine[[-5]]<>"_L_"<>$CommandLine[[-4]]<>"_tau_"<>$CommandLine[[-3]]<>"_N_"<>$CommandLine[[-2]]<>"_kmax_"<>$CommandLine[[-1]]<>"_nat_tau/ClTT_rei_tensor/";
dumpdir="/tomerv3/Sida/b-mode-log-griding/f_"<>$CommandLine[[-5]]<>"_L_"<>$CommandLine[[-4]]<>"_tau_"<>$CommandLine[[-3]]<>"_N_"<>$CommandLine[[-2]]<>"_kmax_"<>$CommandLine[[-1]]<>"_nat_tau/";

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

lmode=Table[i,{i,2,9}]~Join~Table[10i,{i,1,10}];
lmode={200,300,400,500};
lmode={600,700};
lmode={2000};
lmode={700,2000};
lmode={1000};
lmode={500,1000};
lmode={50,80,100,300};
lmode={2,5,8,10,30,2000};
lmode={2,5,8,10,30,50,80,100,300,500,800,1000,2000};

ans=Table[0,Length[lmode]];

allnum=Table[True,Length[lmode]];

For[uu=1,uu<=Ntot(*uu<=93*),uu++,
ans2=Table[0,Length[lmode]];
For[vv=1,vv<=Ntot,vv++,
tab=Import[workdir<>"uu_"<>ToString[uu]<>"_GF/uu_"<>ToString[uu]<>"_vv_"<>ToString[vv]<>"_6.dat","Table"];
If[tab==$Failed,Continue[]];
tab=Flatten[tab];
If[NumberQ[#]&/@tab != allnum,Print["not num: ",uu," ",vv];Continue[]];
ans2+=tab;
ans+=tab;
];
Print[uu," ",ans2];
]
Print[ans]
(* for BM1, 3 for ll=700 2000 *)
(* 5 for ll=1000 *)
Export[dumpdir<>"ClTT_rei_tensor_GF_6.dat",MapThread[{#1,((#1-1)*#1*(#1+1)*(#1+2))#1*(#1+1)*#2/(2 Pi)}&,{lmode,ans}]]
