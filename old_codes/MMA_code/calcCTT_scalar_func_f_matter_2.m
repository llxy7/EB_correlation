Off[NIntegrate::ncvb,NIntegrate::slwcon,NIntegrate::eincr,General::munfl,NIntegrate::izero];

ll=ToExpression[$CommandLine[[-2]]];
kk=ToExpression[$CommandLine[[-1]]];

log10f=14.96;
log10Lambda=ToExpression[$CommandLine[[-6]]];
pushback=ToExpression[$CommandLine[[-5]]];
Ntot=ToExpression[$CommandLine[[-4]]];
coeffk=ToExpression[$CommandLine[[-3]]];

workdir="/tomerv3/Sida/b-mode-log-griding/f_"<>ToString[log10f]<>"_L_"<>$CommandLine[[-6]]<>"_tau_"<>$CommandLine[[-5]]<>"_N_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_nat_tau/parameters/";
dumpdir="/tomerv3/Sida/b-mode-log-griding/f_"<>ToString[log10f]<>"_L_"<>$CommandLine[[-6]]<>"_tau_"<>$CommandLine[[-5]]<>"_N_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_nat_tau/ClTT_rei_scalar_matter/func_f/";

GeVtog=SetPrecision[1.783*10^-24,100];
GeVtom=SetPrecision[1.973*10^-16,100];
G=SetPrecision[6.71*10^-39,100];
Mpl=1/Sqrt[8π G];
H0=SetPrecision[67.8/(299792458/1000)/(3.086*10^22)*GeVtom,100];

Hz[z_]:=SetPrecision[H0*√((1-0.3-0.3/3400)+0.3*(1+z)^3+0.3/3400*(1+z)^4),100];
τz[zz_]:=NIntegrate[1/Hz[z],{z,zz,∞}(*,WorkingPrecision100*)];

{klist,dk,α,Λ,f,m,τosc}=Import[workdir<>"../parameters.mx"];

klist=klist[[2;;-1]];
dk=dk[[2;;-1]];

tabxgrid=(Import[workdir<>"xgrid_HD.dat","Table"]//Flatten);
xgrid=tabxgrid[[1;;-31;;30]]~Join~{tabxgrid[[-1]]};

func$a=Interpolation[ToExpression[Import["/tomerv3/Sida/b-mode-log-griding/tab_a_tau.dat","Table"]]];
func$ap=Interpolation[ToExpression[Import["/tomerv3/Sida/b-mode-log-griding/tab_ap_tau.dat","Table"]]];
rhoD=Interpolation[Table[{τz[z],(3H0^2 0.3*(1+z)^3)/(8π G)},{z,0,2000,0.5}]];

τ0=τz[0];
τr=4*10^40;
τrei=τz[8];
optdep$rei=0.08;
TT=10^40;

func[tdown_]:=Module[{k,phi},
k=klist[[kk]];
phi=Phi/.(NDSolve[{drho'[t]-TT^2*k^2 v[t]==3Phi'[t],v'[t]+TT*func$ap[t*TT]/func$a[t*TT] v[t]==-Phi[t],TT^2*k^2 Phi[t]+3TT func$ap[t*TT]/func$a[t*TT] Phi'[t]+TT^2*3(func$ap[t*TT]/func$a[t*TT])^2 Phi[t]==-TT^2*4π G func$a[t*TT]^2 rhoD[t*TT]*drho[t],Phi[tdown/TT]==1/TT func$a[tdown]/(3func$ap[tdown]),drho[tdown/TT]==0,v[tdown/TT]==0},{Phi,drho,v},{t,tdown/TT,τ0/TT}(*, PrecisionGoal->10, AccuracyGoal->10*)][[1]]);
SphericalBesselJ[ll,klist[[kk]](τ0-tdown)]*TT*phi[tdown/TT]+NIntegrate[SphericalBesselJ[ll,klist[[kk]](τ0-t)]*phi'[t/TT],{t,tdown,τ0}]
]

(*func[tdown_]:=NIntegrate[SphericalBesselJ[ll,klist[[kk]](τ0-t)]Exp[func$h[tdown]-func$h[t]](klist[[kk]]^2func$a[t]^2+3func$ap[t]^2)/(3func$a[t]func$ap[t]),{t,tdown,τ0}(*,WorkingPrecision100*)];*)

Export[dumpdir<>"ll_"<>$CommandLine[[-2]]<>"_kk_"<>$CommandLine[[-1]]<>"_2.dat",Table[{xgrid[[i]],func[xgrid[[i]]]},{i,1,Length[xgrid]}]];
