$MinPrecision=10;
Off[NIntegrate::ncvb,NIntegrate::slwcon]

log10f=ToExpression[$CommandLine[[-5]]];
log10Lambda=ToExpression[$CommandLine[[-4]]];
pushback=ToExpression[$CommandLine[[-3]]];
Ntot=ToExpression[$CommandLine[[-2]]];
coeffk=ToExpression[$CommandLine[[-1]]];

workdir="/tomerv3/Sida/b-mode-log-griding/f_"<>$CommandLine[[-5]]<>"_L_"<>$CommandLine[[-4]]<>"_tau_"<>$CommandLine[[-3]]<>"_N_"<>$CommandLine[[-2]]<>"_kmax_"<>$CommandLine[[-1]]<>"_nat_tau/ndsolve/";
CreateDirectory[workdir]

GeVtog=SetPrecision[1.783*10^-24,100];
GeVtom=SetPrecision[1.973*10^-16,100];
G=SetPrecision[6.71*10^-39,100];
Mpl=1/Sqrt[8π G];
H0=SetPrecision[67.8/(299792458/1000)/(3.086*10^22)*GeVtom,100];
Hz[z_]:=SetPrecision[H0*Sqrt[(1-0.3-0.3/3400)+0.3*(1+z)^3+0.3/3400*(1+z)^4],100];
tz[zz_]:=NIntegrate[1/((1+z)Hz[z]),{z,zz,∞}];(* t as function of z *)
τz[zz_]:=NIntegrate[1/Hz[z],{z,zz,∞},WorkingPrecision->100];
α=400;

τ0=τz[0];
τr=4*10^40;
τosc=5*10^40;
τrei=τz[8];

Λ=10^log10Lambda;
f=500*10^log10f;
m=Λ^2/f/Sqrt[2];

kmax=coeffk*10^(-39);
τmin=τosc;

τr=4*10^40;
tmax=τ0;

(*
explist=Subdivide[Log10[1/(2 τ0)], Log10[kmax], Ntot];
*)
explist=Subdivide[Log10[H0/4], Log10[kmax], Ntot];
klist=(10^#)&/@explist;
Ntot=Ntot+1;
dk={klist[[1]]}~Join~Table[klist[[i]]-klist[[i-1]],{i,2,Ntot}];

aH=Interpolation[Table[{τz[z]//Quiet, 1/(1+z)*Hz[z]},{z,0,2000}]];
a2=Interpolation[Table[{τz[z]//Quiet, (1/(1+z))^2},{z,0,2000}]];

Clear[varlist,eqns,ic];
varlist=Table[Subscript[vr, i,pp],{i,1,Ntot}]~Join~Table[Subscript[vi, i,pp],{i,1,Ntot}]~Join~{ϕ};
eqns=Table[D[Subscript[vr, i,pp][τ],{τ,2}]+((klist[[i]])^2-klist[[i]] α/f ϕ'[τ])Subscript[vr, i,pp][τ]==0,{i,1,Ntot}]~Join~Table[D[Subscript[vi, i,pp][τ],{τ,2}]+((klist[[i]])^2-klist[[i]] α/f ϕ'[τ])Subscript[vi, i,pp][τ]==0,{i,1,Ntot}]~Join~
{ϕ''[τ]+2*aH[τ] ϕ'[τ]+a2[τ] m^2 ϕ[τ]==-(α/f) 1/(2π^2 a2[τ]) Sum[dk[[i]](klist[[i]])^3 (Subscript[vr, i,pp][τ]*Subscript[vr, i,pp]'[τ]+Subscript[vi, i,pp][τ]*Subscript[vi, i,pp]'[τ]),{i,1,Ntot}]};
ic=Flatten[
Table[{Subscript[vr, i,pp][τmin]==Re[Exp[-I klist[[i]] *τmin]/Sqrt[2klist[[i]]]],(D[Subscript[vr, i,pp][x],x]/.{x->τmin})==Re[-I klist[[i]] Exp[-I klist[[i]]*τmin]/Sqrt[2klist[[i]]]]},{i,1,Ntot}]~Join~Table[{Subscript[vi, i,pp][τmin]==Im[Exp[-I klist[[i]] *τmin]/Sqrt[2klist[[i]]]],(D[Subscript[vi, i,pp][x],x]/.{x->τmin})==Im[-I klist[[i]] Exp[-I klist[[i]]*τmin]/Sqrt[2klist[[i]]]]},{i,1,Ntot}]~Join~
{ϕ[τmin]==-f,ϕ'[τmin]==0}];

sol4=NDSolve[eqns~Join~ic,varlist,{τ,τmin,tmax}][[1]];

For[ii=1,ii<=Ntot,Export[workdir<>"vrpp"<>ToString[ii]<>".dat",MapThread[{#1[[1]],#2}&,{(Subscript[vr, ii,pp]/.sol4)["Grid"],(Subscript[vr, ii,pp]/.sol4)["ValuesOnGrid"]}]];
Export[workdir<>"vipp"<>ToString[ii]<>".dat",MapThread[{#1[[1]],#2}&,{(Subscript[vi, ii,pp]/.sol4)["Grid"],(Subscript[vi, ii,pp]/.sol4)["ValuesOnGrid"]}]];
ii++]

Export[workdir<>"phi.dat",MapThread[{#1[[1]],#2}&,{(ϕ/.sol4)["Grid"],(ϕ/.sol4)["ValuesOnGrid"]}]];

(*
Export[workdir<>"../parameters.mx",{klist,dk,Λ,f,m,τosc}];
*)
Export[workdir<>"../parameters.mx",{klist,dk,α,Λ,f,m,τosc}];
