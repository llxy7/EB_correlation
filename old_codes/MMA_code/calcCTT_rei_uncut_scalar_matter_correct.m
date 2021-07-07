$MinPrecision=100;
Off[NIntegrate::ncvb,NIntegrate::slwcon,NIntegrate::eincr,NIntegrate::izero,General::munfl,NIntegrate::inumr];

uu=ToExpression[$CommandLine[[-2]]];
vv=ToExpression[$CommandLine[[-1]]];

log10f=14.96;
log10Lambda=ToExpression[$CommandLine[[-6]]];
pushback=ToExpression[$CommandLine[[-5]]];
Ntot=ToExpression[$CommandLine[[-4]]];
coeffk=ToExpression[$CommandLine[[-3]]];

workdir="/tomerv3/Sida/b-mode-log-griding/f_"<>ToString[log10f]<>"_L_"<>$CommandLine[[-6]]<>"_tau_"<>$CommandLine[[-5]]<>"_N_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_nat_tau/ndsolve/";
dumpdir="/tomerv3/Sida/b-mode-log-griding/f_"<>ToString[log10f]<>"_L_"<>$CommandLine[[-6]]<>"_tau_"<>$CommandLine[[-5]]<>"_N_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_nat_tau/ClTT_rei_scalar_matter/uu_"<>$CommandLine[[-2]]<>"/";

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

func$aH=Interpolation[Import["/tomerv3/Sida/b-mode-log-griding/tab_aH_tau.dat","Table"]];
func$a=Interpolation[Import["/tomerv3/Sida/b-mode-log-griding/tab_a_tau.dat","Table"]];
func$ap=Interpolation[Import["/tomerv3/Sida/b-mode-log-griding/tab_ap_tau.dat","Table"]];

Hz[z_]:=SetPrecision[H0*Sqrt[(1-0.3-0.3/3400)+0.3*(1+z)^3+0.3/3400*(1+z)^4],100];
τz[zz_]:=NIntegrate[1/Hz[z],{z,zz,∞},WorkingPrecision->100];

xgrid=Import[workdir<>"../parameters/xgrid_HD.dat","Table"]//Flatten;
xgrid=xgrid[[1;;-31;;30]]~Join~{xgrid[[-1]]};
tabzero={#,0}&/@xgrid;

τ0=τz[0];
τr=4*10^40;
τrei=τz[8];
optdep$rei=0.08;

vruupp=Interpolation[Import[workdir<>"vrpp"<>ToString[uu+1]<>".dat","Table"]];
viuupp=Interpolation[Import[workdir<>"vipp"<>ToString[uu+1]<>".dat","Table"]];
vrvvpp=Interpolation[Import[workdir<>"vrpp"<>ToString[vv+1]<>".dat","Table"]];
vivvpp=Interpolation[Import[workdir<>"vipp"<>ToString[vv+1]<>".dat","Table"]];
vrwwpp=Interpolation[Import[workdir<>"vrpp"<>ToString[uu+1]<>".dat","Table"]];
viwwpp=Interpolation[Import[workdir<>"vipp"<>ToString[uu+1]<>".dat","Table"]];

lmode={2,5,8,10,30,50,80,100,300,500};

ans=Table[0,Length[lmode]];
wwmin=FirstPosition[klist,SelectFirst[klist,#>=Abs[klist[[uu]]-klist[[vv]]]&]][[1]];
wwmax=FirstPosition[klist,SelectFirst[Reverse[klist],#<=(klist[[uu]]+klist[[vv]])&]][[1]];

tab$func=Table[$Failed,Length[lmode]];
Do[tabtmp=Import[dumpdir<>"../func_f/ll_"<>ToString[lmode[[i]]]<>"_kk_"<>ToString[uu]<>"_2.dat","Table"]//Quiet;
If[tabtmp==$Failed,tabtmp={}];

If[Length[tabtmp]<Length[xgrid],tabtmp=MapThread[{#1,#2}&,{xgrid,(#[[2]]&/@tabtmp)~Join~Table[0,Length[xgrid]-Length[tabtmp]]}]];
tabtmp=DeleteDuplicatesBy[tabtmp,First];

tab$func[[i]]=Interpolation[tabtmp,InterpolationOrder->1]
,{i,1,Length[lmode]}
];

TT=10^40;

Ar[ll_]:=NIntegrate[1/(func$a[t]^2) tab$func[[ll]][t](vrvvpp[t]vrwwpp[t]-vivvpp[t]viwwpp[t]),{t,τosc,τ0}];
Ai[ll_]:=NIntegrate[1/(func$a[t]^2) tab$func[[ll]][t](vrvvpp[t]viwwpp[t]+vivvpp[t]vrwwpp[t]),{t,τosc,τ0}];
Br[ll_]:=NIntegrate[1/(func$a[t]^2) tab$func[[ll]][t](vrvvpp'[t]vrwwpp'[t]-vivvpp'[t]viwwpp'[t]),{t,τosc,τ0}];
Bi[ll_]:=NIntegrate[1/(func$a[t]^2) tab$func[[ll]][t](vrvvpp'[t]viwwpp'[t]+vivvpp'[t]vrwwpp'[t]),{t,τosc,τ0}];

Do[vrwwpp=Interpolation[Import[workdir<>"vrpp"<>ToString[ww+1]<>".dat","Table"]];
viwwpp=Interpolation[Import[workdir<>"vipp"<>ToString[ww+1]<>".dat","Table"]];
tabAr=Ar[#]&/@Table[i,{i,1,Length[lmode]}];
tabAi=Ai[#]&/@Table[i,{i,1,Length[lmode]}];
tabBr=Br[#]&/@Table[i,{i,1,Length[lmode]}];
tabBi=Bi[#]&/@Table[i,{i,1,Length[lmode]}];
ans+=2048π^7 G^2 (klist[[uu]]dk[[uu]]klist[[vv]]dk[[vv]]klist[[ww]]dk[[ww]])/(4π^2) (klist[[vv]]^2+2klist[[vv]]klist[[ww]]+klist[[ww]]^2-klist[[uu]]^2)^2/8*((tabAr+tabBr/(klist[[vv]]klist[[ww]]))^2+(tabAi+tabBi/(klist[[vv]]klist[[ww]]))^2),{ww,wwmin,wwmax}];

(* 2 for ll=200 300 500 *)
(* 3 for ll=800 *)
(* 4 for ll=2000 *)
Export[dumpdir<>"uu_"<>$CommandLine[[-2]]<>"_vv_"<>$CommandLine[[-1]]<>"_correct.dat",ans,"LineSeparators"->" "];
