$MinPrecision=100;
Off[NIntegrate::ncvb,NIntegrate::slwcon,NIntegrate::eincr,General::munfl,NIntegrate::izero];

uu=ToExpression[$CommandLine[[-2]]];
vv=ToExpression[$CommandLine[[-1]]];

log10f=14.96;
(*log10f=14;*)
(*log10f=18;*)
(*log10f=ToExpression[$CommandLine[[-7]]];*)
log10Lambda=ToExpression[$CommandLine[[-6]]];
pushback=ToExpression[$CommandLine[[-5]]];
Ntot=ToExpression[$CommandLine[[-4]]];
coeffk=ToExpression[$CommandLine[[-3]]];

workdir="/tomerv3/Sida/b-mode-log-griding/f_"<>ToString[log10f]<>"_L_"<>$CommandLine[[-6]]<>"_tau_"<>$CommandLine[[-5]]<>"_N_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_nat_tau/ndsolve/";
dumpdir="/tomerv3/Sida/b-mode-log-griding/f_"<>ToString[log10f]<>"_L_"<>$CommandLine[[-6]]<>"_tau_"<>$CommandLine[[-5]]<>"_N_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_nat_tau/ClTT_rei_tensor/uu_"<>$CommandLine[[-2]]<>"_GF/";

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
func$a=Interpolation[ToExpression[Import["/tomerv3/Sida/b-mode-log-griding/tab_a_tau.dat","Table"]]];

Hz[z_]:=SetPrecision[H0*Sqrt[(1-0.3-0.3/3400)+0.3*(1+z)^3+0.3/3400*(1+z)^4],100];
τz[zz_]:=NIntegrate[1/Hz[z],{z,zz,∞},WorkingPrecision->100];

τ0=τz[0];
τr=4*10^40;
τrei=τz[8];
optdep$rei=0.08;

lmode=Table[i,{i,2,10}]~Join~Table[10i,{i,2,10}];
lmode={200,300,400,500};
lmode={600,700};
lmode={2000};
lmode={700,2000};
lmode={1000};
lmode={500,1000};
lmode={50,80,100,300};
lmode={2,5,8,10,30,2000};
lmode={2,5,8,10,30,50,80,100,300,500,800,1000,2000};

xgridfull=Import[workdir<>"../parameters/xgrid_HD.dat","Table"];
xgrid=xgridfull[[1;;-11;;10]]~Join~{xgridfull[[-1]]};

func=Table[$Failed,Length[lmode]];
For[i=1,i<=Length[lmode],i++,
tabfunc=Import[workdir<>"../ClTT_rei_tensor/func_f_GF/ll_"<>ToString[lmode[[i]]]<>"_kk_"<>ToString[uu]<>"_2.dat","Table"]//Quiet;
If[tabfunc==$Failed,Continue[]];
(*If[Length[tabfunc]<Length[xgrid],tabfunc=MapThread[{#1,#2}&,{xgrid,(#[[2]]&/@tabfunc)~Join~Table[0,Length[xgrid]-Length[tabfunc]]}]];*)
tabfunc=DeleteDuplicatesBy[tabfunc,First];
func[[i]]=Interpolation[tabfunc,InterpolationOrder->1]
]

vruupp=Interpolation[Import[workdir<>"vrpp"<>ToString[uu+1]<>".dat","Table"]];
viuupp=Interpolation[Import[workdir<>"vipp"<>ToString[uu+1]<>".dat","Table"]];
vrvvpp=Interpolation[Import[workdir<>"vrpp"<>ToString[vv+1]<>".dat","Table"]];
vivvpp=Interpolation[Import[workdir<>"vipp"<>ToString[vv+1]<>".dat","Table"]];
vrwwpp=Interpolation[Import[workdir<>"vrpp"<>ToString[uu+1]<>".dat","Table"]];
viwwpp=Interpolation[Import[workdir<>"vipp"<>ToString[uu+1]<>".dat","Table"]];

Θ[a_,b_,c_]=1/16 ((1+(a^2+b^2-c^2)/(2a b))^2 (1+(a^2-b^2+c^2)/(2a c))^2+(1-(a^2+b^2-c^2)/(2a b))^2 (1-(a^2-b^2+c^2)/(2a c))^2);

Vr[l$index_,ww_]:=NIntegrate[2/Mpl^2 1/func$a[t] (func[[l$index]][t])(klist[[vv]]klist[[ww]](vrvvpp[t]vrwwpp[t]-vivvpp[t]viwwpp[t])+(vrvvpp'[t]vrwwpp'[t]-vivvpp'[t]viwwpp'[t])),{t,τosc,999τ0/1000},WorkingPrecision->100]//Re
Vi[l$index_,ww_]:=NIntegrate[2/Mpl^2 1/func$a[t] (func[[l$index]][t])(klist[[vv]]klist[[ww]](vrvvpp[t]viwwpp[t]+vivvpp[t]vrwwpp[t])+(vrvvpp'[t]viwwpp'[t]+vivvpp'[t]vrwwpp'[t])),{t,τosc,999τ0/1000},WorkingPrecision->100]//Re

ans=Table[0,Length[lmode]];
wwmin=FirstPosition[klist,SelectFirst[klist,#>=Abs[klist[[uu]]-klist[[vv]]]&]][[1]];
wwmax=FirstPosition[klist,SelectFirst[Reverse[klist],#<=(klist[[uu]]+klist[[vv]])&]][[1]];

calcClTT[l$index_,ww_]:=1/(4π^5) (9π)/2(*36π^2optdep$rei^2*)klist[[uu]]dk[[uu]]klist[[vv]]dk[[vv]]klist[[ww]]dk[[ww]]Θ[klist[[uu]],klist[[vv]],klist[[ww]]](Vr[l$index,ww]^2+Vi[l$index,ww]^2);

Do[vrwwpp=Interpolation[Import[workdir<>"vrpp"<>ToString[ww+1]<>".dat","Table"]];
viwwpp=Interpolation[Import[workdir<>"vipp"<>ToString[ww+1]<>".dat","Table"]];
For[lind=1,lind<=Length[lmode],lind++,
If[func[[lind]]==$Failed,Continue[]];
ans[[lind]]+=calcClTT[lind,ww]
]
,{ww,wwmin,wwmax}];

(* 3 for l=600, 700 *)
(* 4 for l=2000 *)
(* for BM1, 3 for ll=700 2000 *)
(* 5 for l=1000 *)
Export[dumpdir<>"uu_"<>$CommandLine[[-2]]<>"_vv_"<>$CommandLine[[-1]]<>"_6.dat",ans,"LineSeparators"->" "];

