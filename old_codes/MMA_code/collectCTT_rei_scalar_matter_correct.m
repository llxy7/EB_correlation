$MinPrecision=100;
Off[NIntegrate::ncvb,NIntegrate::slwcon,NIntegrate::eincr,General::munfl];

log10f=ToExpression[$CommandLine[[-5]]];
log10Lambda=ToExpression[$CommandLine[[-4]]];
pushback=ToExpression[$CommandLine[[-3]]];
Ntot=ToExpression[$CommandLine[[-2]]];
coeffk=ToExpression[$CommandLine[[-1]]];

workdir="/tomerv3/Sida/b-mode-log-griding/f_"<>$CommandLine[[-5]]<>"_L_"<>$CommandLine[[-4]]<>"_tau_"<>$CommandLine[[-3]]<>"_N_"<>$CommandLine[[-2]]<>"_kmax_"<>$CommandLine[[-1]]<>"_nat_tau/ClTT_rei_scalar_matter/";
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
lmode={2,3,5,8,10,20,30,50,70,100};
lmode={200,300,500};
lmode={800};
(*
lmode={200,300,400,500};
lmode={70,80,90,100};
lmode={600,700};
lmode={2000};
*)

ans=Table[0,Length[lmode]];
ans2=Table[0,Length[lmode]];
ans3=Table[0,Length[lmode]];

allnum=Table[True,Length[lmode]];

For[uu=1,uu<=Ntot,uu++,
ans2=Table[0,Length[lmode]];
For[vv=1,vv<=Ntot,vv++,
impt1=Import[workdir<>"uu_"<>ToString[uu]<>"/uu_"<>ToString[uu]<>"_vv_"<>ToString[vv]<>"_correct_3.dat","Table"];
If[impt1==$Failed, Print[ToString[uu]<>" "<>ToString[vv]];Continue[]];
impt1=Flatten[impt1];
If[(NumberQ[#]&/@impt1)!=allnum,Print["not number "<>ToString[uu]<>" "<>ToString[vv]];Continue[]];
ans+=Flatten[impt1];
ans2+=Flatten[impt1];
];
Print[uu," ",ans2];
]
Print[ans]
Print[MapThread[{#1,#1*(#1+1)*#2/(2 Pi)}&,{lmode,ans}]]

TT=10^40;

(* 2 for l=200 300 500 *)
(* 3 for l=800 *)
(* 5 for l=2000 *)
Export[dumpdir<>"ClTT_rei_scalar_matter_correct_3.dat",MapThread[{#1,#1*(#1+1)*#2/(2 Pi)}&,{lmode,ans}]]

