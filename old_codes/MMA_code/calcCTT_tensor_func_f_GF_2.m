Off[NIntegrate::ncvb,NIntegrate::slwcon];

ll=ToExpression[$CommandLine[[-2]]];
kk=ToExpression[$CommandLine[[-1]]];

log10f=14.96;
log10Lambda=ToExpression[$CommandLine[[-6]]];
pushback=ToExpression[$CommandLine[[-5]]];
Ntot=ToExpression[$CommandLine[[-4]]];
coeffk=ToExpression[$CommandLine[[-3]]];

workdir="/tomerv3/Sida/b-mode-log-griding/f_"<>ToString[log10f]<>"_L_"<>$CommandLine[[-6]]<>"_tau_"<>$CommandLine[[-5]]<>"_N_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_nat_tau/parameters/";
dumpdir="/tomerv3/Sida/b-mode-log-griding/f_"<>ToString[log10f]<>"_L_"<>$CommandLine[[-6]]<>"_tau_"<>$CommandLine[[-5]]<>"_N_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_nat_tau/ClTT_rei_tensor/func_f_GF/";

GeVtog=SetPrecision[1.783*10^-24,100];
GeVtom=SetPrecision[1.973*10^-16,100];
G=SetPrecision[6.71*10^-39,100];
Mpl=1/Sqrt[8π G];
H0=SetPrecision[67.8/(299792458/1000)/(3.086*10^22)*GeVtom,100];

Hz[z_]:=SetPrecision[H0*√((1-0.3-0.3/3400)+0.3*(1+z)^3+0.3/3400*(1+z)^4),100];
τz[zz_]:=NIntegrate[1/Hz[z],{z,zz,∞},WorkingPrecision->100];

{klist,dk,α,Λ,f,m,τosc}=Import[workdir<>"../parameters.mx"];

klist=klist[[2;;-1]];
dk=dk[[2;;-1]];

klist=SetPrecision[klist,100];
dk=SetPrecision[dk,100];

τ0=τz[0];
τr=4*10^40;
τrei=τz[8];
optdep$rei=0.08;

tabxgrid=(Import[workdir<>"xgrid_HD.dat","Table"]//Flatten);
(*xgrid=tabxgrid[[1;;-11;;10]];*)
xgrid=tabxgrid[[1;;-31;;30]];

tabatau=Import["/tomerv3/Sida/b-mode-log-griding/tab_a_tau.dat","Table"];
tabaptau=Import["/tomerv3/Sida/b-mode-log-griding/tab_ap_tau.dat","Table"];
afit=Fit[tabatau,{1,t,t^2,t^3,t^4,t^5,t^6},t];
apfit=Fit[tabaptau,{1,t,t^2,t^3,t^4,t^5,t^6},t];

func$a=Interpolation[tabatau];
func$ap=Interpolation[tabaptau];

func$a$num[t_?NumberQ]:=func$a[t];
func$ap$num[t_?NumberQ]:=func$ap[t];

g1=Interpolation[Import[workdir<>"../GF_new/g1_kk_"<>$CommandLine[[-1]]<>".dat"]];
g1p=Interpolation[Import[workdir<>"../GF_new/g1p_kk_"<>$CommandLine[[-1]]<>".dat"]];
g2=Interpolation[Import[workdir<>"../GF_new/g2_kk_"<>$CommandLine[[-1]]<>".dat"]];
g2p=Interpolation[Import[workdir<>"../GF_new/g2p_kk_"<>$CommandLine[[-1]]<>".dat"]];

g1$num[t_?NumberQ]:=g1[t];
g1p$num[t_?NumberQ]:=g1p[t];
g2$num[t_?NumberQ]:=g2[t];
g2p$num[t_?NumberQ]:=g2p[t];

(* G=g2[tp]/(g1p[tp]g2[tp]-g1[tp]g2p[tp])g1[t]-g1[tp]/(g1p[tp]g2[tp]-g1[tp]g2p[tp])g2[t] *)
func[tp_]:=NIntegrate[1/func$a$num[t] (g2$num[tp]g1p$num[t]-g1$num[tp]g2p$num[t]-func$ap$num[t]/func$a$num[t] (g2$num[tp]g1$num[t]-g1$num[tp]g2$num[t]))SphericalBesselJ[ll,(τ0-t)klist[[kk]]] 1/((τ0-t)^2 klist[[kk]]^2),{t,tp,τ0}(*,WorkingPrecision100*)];
(*
For[i=1,i≤Length[xgrid],i++,
Print[i," ",xgrid[[i]]," ",func[xgrid[[i]]]]
]
*)

ans=({#,func[#]}&/@xgrid)~Join~{{N[τ0,10],0}};

Export[dumpdir<>"ll_"<>$CommandLine[[-2]]<>"_kk_"<>$CommandLine[[-1]]<>"_2.dat",ans];
