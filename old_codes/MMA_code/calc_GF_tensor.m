kk=ToExpression[$CommandLine[[-1]]];

log10f=ToExpression[$CommandLine[[-6]]];;
log10Lambda=ToExpression[$CommandLine[[-5]]];
pushback=ToExpression[$CommandLine[[-4]]];
Ntot=ToExpression[$CommandLine[[-3]]];
coeffk=ToExpression[$CommandLine[[-2]]];

workdir="/tomerv3/Sida/b-mode-log-griding/f_"<>$CommandLine[[-6]]<>"_L_"<>$CommandLine[[-5]]<>"_tau_"<>$CommandLine[[-4]]<>"_N_"<>$CommandLine[[-3]]<>"_kmax_"<>$CommandLine[[-2]]<>"_nat_tau/parameters/";
dumpdir="/tomerv3/Sida/b-mode-log-griding/f_"<>$CommandLine[[-6]]<>"_L_"<>$CommandLine[[-5]]<>"_tau_"<>$CommandLine[[-4]]<>"_N_"<>$CommandLine[[-3]]<>"_kmax_"<>$CommandLine[[-2]]<>"_nat_tau/GF_new/";
If[DirectoryQ[dumpdir]!=True,CreateDirectory[dumpdir]];

GeVtog=SetPrecision[1.783*10^-24,100];
GeVtom=SetPrecision[1.973*10^-16,100];
G=SetPrecision[6.71*10^-39,100];
Mpl=1/Sqrt[8π G];
H0=SetPrecision[67.8/(299792458/1000)/(3.086*10^22)*GeVtom,100];

{klist,dk,α,Λ,f,m,τosc}=Import[workdir<>"../parameters.mx"];

klist=klist[[2;;-1]];
dk=dk[[2;;-1]];

klist=SetPrecision[klist,100];
dk=SetPrecision[dk,100];

func$appoa=Interpolation[Import["/tomerv3/Sida/b-mode-log-griding/tab_appoa_tau.dat","Table"]];
func$appoa$num[t_?NumberQ]:=func$appoa[t];

Hz[z_]:=SetPrecision[H0*Sqrt[(1-0.3-0.3/3400)+0.3*(1+z)^3+0.3/3400*(1+z)^4],100];
τz[zz_]:=NIntegrate[1/Hz[z],{z,zz,∞},WorkingPrecision->100];

τ0=τz[0];
τr=4*10^40;
τrei=τz[8];
optdep$rei=0.08;

xgrid=(Import[workdir<>"xgrid_HD.dat","Table"]//Flatten);

x0=10^40;

g1=g/.NDSolve[{g''[t]+(klist[[kk]]^2-func$appoa$num[t])g[t]==0,g[τosc]==0,g'[τosc]==1/x0},g,{t,τosc,τ0},MaxSteps->10^8][[1]]
g2=g/.NDSolve[{g''[t]+(klist[[kk]]^2-func$appoa$num[t])g[t]==0,g[τosc]==1,g'[τosc]==0},g,{t,τosc,τ0},MaxSteps->10^8][[1]]

Export[dumpdir<>"g1_kk_"<>$CommandLine[[-1]]<>".dat",{#,x0*g1[#]}&/@xgrid];
Export[dumpdir<>"g1p_kk_"<>$CommandLine[[-1]]<>".dat",{#,x0*g1'[#]}&/@xgrid];
Export[dumpdir<>"g2_kk_"<>$CommandLine[[-1]]<>".dat",{#,g2[#]}&/@xgrid];
Export[dumpdir<>"g2p_kk_"<>$CommandLine[[-1]]<>".dat",{#,g2'[#]}&/@xgrid];
