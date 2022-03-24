$MinPrecision=10;
Off[NIntegrate::ncvb,NIntegrate::slwcon,NDSolve::precw]

wp=30;
log10f=SetPrecision[ToExpression[$CommandLine[[-4]]],wp];
log10Lambda=SetPrecision[ToExpression[$CommandLine[[-3]]],wp];
Ntot=200;
coeffk=ToExpression[$CommandLine[[-2]]];
α=ToExpression[$CommandLine[[-1]]];

workdir="/tomerv3/Sida/EB_correlation/f_"<>$CommandLine[[-4]]<>"_L_"<>$CommandLine[[-3]]<>"_kmax_"<>$CommandLine[[-2]]<>"_alpha_"<>$CommandLine[[-1]]<>"/ndsolve/";
CreateDirectory[workdir]

GeVtog=SetPrecision[1.783*10^-24,wp];
GeVtom=SetPrecision[1.973*10^-16,wp];
GeVtoMpc=SetPrecision[GeVtom/(3.086*10^22),wp];
G=SetPrecision[6.71*10^-39*GeVtoMpc^2,wp];
Mpl=1/Sqrt[8π G];
H0=SetPrecision[67.8/(299792458/1000)/(3.086*10^22)*GeVtom/GeVtoMpc,wp];
Hz[z_]:=SetPrecision[H0*Sqrt[(1-0.3-0.3/3400)+0.3*(1+z)^3+0.3/3400*(1+z)^4],wp];
τz[zz_]:=NIntegrate[1/Hz[z],{z,zz,∞},WorkingPrecision->wp];

(*α=400;*)

τ0=τz[0];
τosc=5*10^40*GeVtoMpc;
τosc=200;
τrei=τz[8];

Λ=SetPrecision[10^log10Lambda/GeVtoMpc,wp];
f=SetPrecision[500*10^log10f/GeVtoMpc,wp];
m=Λ^2/f/Sqrt[2];

kmax=coeffk*10^(-39)/GeVtoMpc;
τmin=τosc;

tmax=τ0;

explist=Subdivide[Log10[H0/4], Log10[kmax], Ntot];
klist=(10^#)&/@explist;
Ntot=Ntot+1;
dk={klist[[1]]}~Join~Table[klist[[i]]-klist[[i-1]],{i,2,Ntot}];

aH=Interpolation[Import["/tomerv3/Sida/EB_correlation/tab_aH_tau.dat"]];
a2=Interpolation[{#[[1]],#[[2]]^2}&/@Import["/tomerv3/Sida/EB_correlation/tab_a_tau.dat"]];
(*aH=Interpolation[Table[{τz[z]//Quiet, 1/(1+z)*Hz[z]},{z,0,2000,1/2}]];
a2=Interpolation[Table[{τz[z]//Quiet, (1/(1+z))^2},{z,0,2000,1/2}]];*)

Clear[varlist,eqns,ic];
varlist=Table[Subscript[vr, i,pp],{i,1,Ntot}]~Join~Table[Subscript[vi, i,pp],{i,1,Ntot}]~Join~{ϕ};
eqns=Table[D[Subscript[vr, i,pp][τ],{τ,2}]+((klist[[i]])^2-klist[[i]] α ϕ'[τ])Subscript[vr, i,pp][τ]==0,{i,1,Ntot}]~Join~Table[D[Subscript[vi, i,pp][τ],{τ,2}]+((klist[[i]])^2-klist[[i]] α ϕ'[τ])Subscript[vi, i,pp][τ]==0,{i,1,Ntot}]~Join~
{ϕ''[τ]+2*aH[τ] ϕ'[τ]+a2[τ] m^2 ϕ[τ]==-(α/f^2) 1/(2π^2 a2[τ]) Sum[dk[[i]](klist[[i]])^3 (Subscript[vr, i,pp][τ]*Subscript[vr, i,pp]'[τ]+Subscript[vi, i,pp][τ]*Subscript[vi, i,pp]'[τ]),{i,1,Ntot}]};
ic=Flatten[
Table[{Subscript[vr, i,pp][τmin]==Re[Exp[-I klist[[i]] *τmin]/Sqrt[2klist[[i]]]],(D[Subscript[vr, i,pp][x],x]/.{x->τmin})==Re[-I klist[[i]] Exp[-I klist[[i]]*τmin]/Sqrt[2klist[[i]]]]},{i,1,Ntot}]~Join~Table[{Subscript[vi, i,pp][τmin]==Im[Exp[-I klist[[i]] *τmin]/Sqrt[2klist[[i]]]],(D[Subscript[vi, i,pp][x],x]/.{x->τmin})==Im[-I klist[[i]] Exp[-I klist[[i]]*τmin]/Sqrt[2klist[[i]]]]},{i,1,Ntot}]~Join~
{ϕ[τmin]==-1,ϕ'[τmin]==0}];

sol4=NDSolve[eqns~Join~ic,varlist,{τ,τmin,tmax}(*,WorkingPrecision->wp*),MaxSteps->10^8][[1]];

Print["solved"];

tabxgrid=Flatten[(ϕ/.sol4)["Grid"]];
Print[Length[tabxgrid]];
tabphi=(ϕ/.sol4)["ValuesOnGrid"];
Export[workdir<>"phi.dat",MapThread[{#1,#2}&, {tabxgrid,tabphi}]];
tabvr=Table[(Subscript[vr, ii,pp]/.sol4)["ValuesOnGrid"],{ii,2,Ntot}];
tabvi=Table[(Subscript[vi, ii,pp]/.sol4)["ValuesOnGrid"],{ii,2,Ntot}];
tabvrp=Table[(Subscript[vr, ii,pp]/.sol4)'["ValuesOnGrid"],{ii,2,Ntot}];
tabvip=Table[(Subscript[vi, ii,pp]/.sol4)'["ValuesOnGrid"],{ii,2,Ntot}];

If[Mod[Length[tabxgrid],10]==1,
tabxgrid=tabxgrid[[1;;-1;;10]];
tabvr=tabvr[[;;, 1;;-1;;10]];
tabvi=tabvi[[;;, 1;;-1;;10]];
tabvrp=tabvrp[[;;, 1;;-1;;10]];
tabvip=tabvip[[;;, 1;;-1;;10]];
,
tabxgrid=tabxgrid[[1;;-1;;10]]~Join~{tabxgrid[[-1]]};
tabvr=Join[tabvr[[;;, 1;;-1;;10]],tabvr[[;;, -1;;-1]],2];
tabvi=Join[tabvi[[;;, 1;;-1;;10]],tabvi[[;;, -1;;-1]],2];
tabvrp=Join[tabvrp[[;;, 1;;-1;;10]],tabvrp[[;;, -1;;-1]],2];
tabvip=Join[tabvip[[;;, 1;;-1;;10]],tabvip[[;;, -1;;-1]],2];
];
DumpSave[workdir<>"sol.mx",{tabxgrid,tabphi,tabvr,tabvi,tabvrp,tabvip}];

Export[workdir<>"../parameters.mx",{klist,dk,α,Λ,f,m,τosc}];
