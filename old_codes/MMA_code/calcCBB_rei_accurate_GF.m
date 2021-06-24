$MinPrecision=100;
Off[NIntegrate::ncvb,NIntegrate::slwcon,NIntegrate::eincr,General::munfl];

uu=ToExpression[$CommandLine[[-2]]];
vv=ToExpression[$CommandLine[[-1]]];

log10f=14.96;
log10Lambda=ToExpression[$CommandLine[[-6]]];
pushback=ToExpression[$CommandLine[[-5]]];
Ntot=ToExpression[$CommandLine[[-4]]];
coeffk=ToExpression[$CommandLine[[-3]]];

workdir="/tomerv3/Sida/b-mode-log-griding/f_"<>ToString[log10f]<>"_L_"<>$CommandLine[[-6]]<>"_tau_"<>$CommandLine[[-5]]<>"_N_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_nat_tau/ndsolve/";
dumpdir="/tomerv3/Sida/b-mode-log-griding/f_"<>ToString[log10f]<>"_L_"<>$CommandLine[[-6]]<>"_tau_"<>$CommandLine[[-5]]<>"_N_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_nat_tau/ClBB_rei_accurate/uu_"<>$CommandLine[[-2]]<>"_GF/";
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

func$aH=Interpolation[Import["/tomerv3/Sida/b-mode-log-griding/tab_aH_tau.dat","Table"]];
func$a=Interpolation[ToExpression[Import["/tomerv3/Sida/b-mode-log-griding/tab_a_tau.dat","Table"]]];

Hz[z_]:=SetPrecision[H0*Sqrt[(1-0.3-0.3/3400)+0.3*(1+z)^3+0.3/3400*(1+z)^4],100];
τz[zz_]:=NIntegrate[1/Hz[z],{z,zz,∞},WorkingPrecision->100];

τ0=τz[0];
τr=4*10^40;
τrei=τz[8];
optdep$rei=0.08;

(*
Print[func$a[τrei]];
*)

vruupp=Interpolation[Import[workdir<>"vrpp"<>ToString[uu+1]<>".dat","Table"]];
viuupp=Interpolation[Import[workdir<>"vipp"<>ToString[uu+1]<>".dat","Table"]];
vrvvpp=Interpolation[Import[workdir<>"vrpp"<>ToString[vv+1]<>".dat","Table"]];
vivvpp=Interpolation[Import[workdir<>"vipp"<>ToString[vv+1]<>".dat","Table"]];
vrwwpp=Interpolation[Import[workdir<>"vrpp"<>ToString[uu+1]<>".dat","Table"]];
viwwpp=Interpolation[Import[workdir<>"vipp"<>ToString[uu+1]<>".dat","Table"]];

Θ[a_,b_,c_]=1/16 ((1+(a^2+b^2-c^2)/(2a b))^2 (1+(a^2-b^2+c^2)/(2a c))^2+(1-(a^2+b^2-c^2)/(2a b))^2 (1-(a^2-b^2+c^2)/(2a c))^2);

funcf=Interpolation[Import[workdir<>"../ClBB_rei_accurate/func_f_GF/kk_"<>$CommandLine[[-2]]<>".dat","Table"],InterpolationOrder->1];

Vr[ww_]:=NIntegrate[2/Mpl^2 1/func$a[t] funcf[t](klist[[vv]]klist[[ww]](vrvvpp[t]vrwwpp[t]-vivvpp[t]viwwpp[t])+(vrvvpp'[t]vrwwpp'[t]-vivvpp'[t]viwwpp'[t])),{t,τosc,999/1000 τrei}, WorkingPrecision->100]//Re
Vi[ww_]:=NIntegrate[2/Mpl^2 1/func$a[t] funcf[t](klist[[vv]]klist[[ww]](vrvvpp[t]viwwpp[t]+vivvpp[t]vrwwpp[t])+(vrvvpp'[t]viwwpp'[t]+vivvpp'[t]vrwwpp'[t])),{t,τosc,999/1000 τrei}, WorkingPrecision->100]//Re

(*
Print[func[H0,klist[[uu]],τrei/2,τrei]];
Print[Vr[uu]];
*)
ans=Table[0,23];
lmode=Table[i,{i,1,9}]~Join~Table[10i,{i,1,9}]~Join~Table[100i,{i,1,5}];
wwmin=FirstPosition[klist,SelectFirst[klist,#>=Abs[klist[[uu]]-klist[[vv]]]&]][[1]];
wwmax=FirstPosition[klist,SelectFirst[Reverse[klist],#<=(klist[[uu]]+klist[[vv]])&]][[1]];

calcClBB[ww_]:=9/π^3 optdep$rei^2 klist[[uu]]dk[[uu]]klist[[vv]]dk[[vv]]klist[[ww]]dk[[ww]]Θ[klist[[uu]],klist[[vv]],klist[[ww]]](Vr[ww]^2+Vi[ww]^2);
tab2=Table[((l+2)/(2l+1) SphericalBesselJ[l-1,klist[[uu]] (τ0-τrei)]-(l-1)/(2l+1) SphericalBesselJ[l+1,klist[[uu]] (τ0-τrei)])^2/.l->lmode[[i]],{i,1,Length[lmode]}];

Do[vrwwpp=Interpolation[Import[workdir<>"vrpp"<>ToString[ww+1]<>".dat","Table"]];
viwwpp=Interpolation[Import[workdir<>"vipp"<>ToString[ww+1]<>".dat","Table"]];
ans+=calcClBB[ww]*tab2,{ww,wwmin,wwmax}];

Export[dumpdir<>"uu_"<>$CommandLine[[-2]]<>"_vv_"<>$CommandLine[[-1]]<>".dat",ans,"LineSeparators"->" "];
