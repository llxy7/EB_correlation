$MinPrecision=10;
Off[NIntegrate::ncvb,NIntegrate::slwcon,General::munfl]

wp=40;
log10f=ToExpression[$CommandLine[[-5]]];
log10Lambda=ToExpression[$CommandLine[[-4]]];
Ntot=200;
coeffk=ToExpression[$CommandLine[[-3]]];
(*alpha=ToExpression[$CommandLine[[-2]]];*)
kk=ToExpression[$CommandLine[[-1]]];
uu=kk;

(* ---------- initialization ---------- *)

workdir="/tomerv3/Sida/EB_correlation/f_"<>$CommandLine[[-5]]<>"_L_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_alpha_"<>$CommandLine[[-2]]<>"/ndsolve/";
dumpdir="/tomerv3/Sida/EB_correlation/f_"<>$CommandLine[[-5]]<>"_L_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_alpha_"<>$CommandLine[[-2]]<>"/ClTT_tensor/";
 
(*{tabxgrid,tabphi,tabvr,tabvi,tabvrp,tabvip}=*)Import[workdir<>"sol.mx"];
{klist,dk,α,Λ,f,m,τosc}=Import[workdir<>"../parameters.mx"];
kmax=coeffk*10^(-39)/GeVtoMpc;

tabxgrid=SetPrecision[tabxgrid,wp];

klist=klist[[2;;-1]];
dk=dk[[2;;-1]];
totpoint=Length[tabxgrid];
intervals=Table[{tabxgrid[[i-1]],tabxgrid[[i]]},{i,2,totpoint}];
intervalLengths=(#[[2]]-#[[1]])&/@intervals;

GeVtog=SetPrecision[1.783*10^-24,wp];
GeVtom=SetPrecision[1.973*10^-16,wp];
GeVtoMpc=SetPrecision[GeVtom/(3.086*10^22),wp];
G=SetPrecision[6.71*10^-39*GeVtoMpc^2,wp];
Mpl=1/Sqrt[8π G];
H0=SetPrecision[67.8/(299792458/1000)/(3.086*10^22)*GeVtom/GeVtoMpc,wp];
Hz[z_]:=SetPrecision[H0*Sqrt[(1-0.3-0.3/3400)+0.3*(1+z)^3+0.3/3400*(1+z)^4],wp];
τz[zz_]:=NIntegrate[1/Hz[z],{z,zz,∞},WorkingPrecision->wp];
optdep$rei=0.08;

τ0=τz[0];
τr=τz[1100];
τrei=τz[8];

τmin=τosc;
tmax=τ0;
Λ=SetPrecision[10^log10Lambda/GeVtoMpc,wp];
f=SetPrecision[500*10^log10f/GeVtoMpc,wp];
m=Λ^2/f/Sqrt[2];

aH=Interpolation[Import["/tomerv3/Sida/EB_correlation/tab_aH_tau.dat","Table"]];
aH$num[t_?NumberQ]:=aH[t];

a=Interpolation[Import["/tomerv3/Sida/EB_correlation/tab_a_tau.dat","Table"]];
a$num[t_?NumberQ]:=a[t];

ap=Interpolation[Import["/tomerv3/Sida/EB_correlation/tab_ap_tau.dat","Table"]];
ap$num[t_?NumberQ]:=ap[t];

appoa=Interpolation[Import["/tomerv3/Sida/EB_correlation/tab_appoa_tau.dat","Table"]];
appoa$num[t_?NumberQ]:=appoa[t];

arecp=Interpolation[Import["/tomerv3/Sida/EB_correlation/tab_1oa_tau.dat","Table"]];
arecp$num[t_?NumberQ]:=arecp[t];

(* ---------- Solving for the GF ---------- *)

g1=g/.(NDSolve[{g''[t]+(klist[[kk]]^2-appoa$num[t])g[t]==0,g[τosc]==0,g'[τosc]==1},g,{t,τosc,τ0},MaxSteps->10^8][[1]]);
g2=g/.(NDSolve[{g''[t]+(klist[[kk]]^2-appoa$num[t])g[t]==0,g[τosc]==1,g'[τosc]==0},g,{t,τosc,τ0},MaxSteps->10^8][[1]]);
g1$num[t_?NumberQ]:=g1[t];
g1p$num[t_?NumberQ]:=g1'[t];
g2$num[t_?NumberQ]:=g2[t];
g2p$num[t_?NumberQ]:=g2'[t];
y$g1=g1$num[#]&/@tabxgrid;
y$g2=g2$num[#]&/@tabxgrid;
y$g1p=g1p$num[#]&/@tabxgrid;
y$g2p=g2p$num[#]&/@tabxgrid;


y$aH=aH$num[#]&/@tabxgrid;
y$arecp=arecp$num[#]&/@tabxgrid;

(* ---------- Tabularizing Spherical Bessel functions, for T ---------- *)

lmode={5,10,50,100,200,300,500,800,2000};

ans=Table[0.,Length[lmode]];

Θ[a_,b_,c_]:=1/16 ((1+(a^2+b^2-c^2)/(2a b))^2 (1+(a^2-b^2+c^2)/(2a c))^2+(1-(a^2+b^2-c^2)/(2a b))^2 (1-(a^2-b^2+c^2)/(2a c))^2);

y$jl=((SphericalBesselJ[lmode,klist[[kk]](τ0-#)]/(klist[[kk]](τ0-#))^2)&/@tabxgrid[[1;;-2]])~Join~{Table[0,Length[lmode]]};

y1=y$arecp*(y$g1p-y$aH*y$g1)*y$jl;
y2=y$arecp*(y$g2p-y$aH*y$g2)*y$jl;
val1=((y1[[1;;-2]]+y1[[2;;-1]])*intervalLengths/2)//Reverse//Accumulate//Reverse;
val2=((y2[[1;;-2]]+y2[[2;;-1]])*intervalLengths/2)//Reverse//Accumulate//Reverse;

fT$tensor=(y$g2[[1;;-2]]*val1-y$g1[[1;;-2]]*val2)~Join~{Table[0,Length[lmode]]};

Do[
wwmin=FirstPosition[klist,SelectFirst[klist,#>=Abs[klist[[uu]]-klist[[vv]]]&]][[1]];wwmax=FirstPosition[klist,SelectFirst[Reverse[klist],#<=(klist[[uu]]+klist[[vv]])&]][[1]];
Do[
intvr=(2/Mpl^2)*y$arecp*fT$tensor*(klist[[vv]]klist[[ww]](tabvr[[vv]]*tabvr[[ww]]-tabvi[[vv]]*tabvi[[ww]])+(tabvrp[[vv]]*tabvrp[[ww]]-tabvip[[vv]]*tabvip[[ww]]));
intvi=(2/Mpl^2)*y$arecp*fT$tensor*(klist[[vv]]klist[[ww]](tabvr[[vv]]*tabvi[[ww]]+tabvi[[vv]]*tabvr[[ww]])+(tabvrp[[vv]]*tabvip[[ww]]+tabvip[[vv]]*tabvrp[[ww]]));
Vr=Total[(intvr[[1;;-2]]+intvr[[2;;-1]])*intervalLengths/2];
Vi=Total[(intvi[[1;;-2]]+intvi[[2;;-1]])*intervalLengths/2];
ans+=1/(4π^5)*(9π/2)*klist[[uu]]dk[[uu]]klist[[vv]]dk[[vv]]klist[[ww]]dk[[ww]]*Θ[klist[[uu]],klist[[vv]],klist[[ww]]](Vr^2+Vi^2);
,{ww,wwmin,wwmax}];
,{vv,1,Ntot}];

Export[dumpdir<>"uu_"<>$CommandLine[[-1]]<>".dat",MapThread[{#1,#2}&,{lmode,ans}]];
(* ClBB *)

(*
dumpdir="/tomerv3/Sida/EB_correlation/f_"<>$CommandLine[[-5]]<>"_L_"<>$CommandLine[[-4]]<>"_kmax_"<>$CommandLine[[-3]]<>"_alpha_"<>$CommandLine[[-2]]<>"/ClBB/";
 
tabxgrid=Select[tabxgrid,#<=τrei&]~Join~{τrei};
totpoint=Length[tabxgrid];
intervals=Table[{tabxgrid[[i-1]],tabxgrid[[i]]},{i,2,totpoint}];
intervalLengths=(#[[2]]-#[[1]])&/@intervals;

y$g1=g1$num[#]&/@tabxgrid;
y$g2=g2$num[#]&/@tabxgrid;
y$g1p=g1p$num[#]&/@tabxgrid;
y$g2p=g2p$num[#]&/@tabxgrid;
y$aH=aH$num[#]&/@tabxgrid;
y$arecp=arecp$num[#]&/@tabxgrid;

y$jl=((SphericalBesselJ[2,klist[[kk]](τrei-#)]/(klist[[kk]](τrei-#))^2)&/@tabxgrid[[1;;-2]])~Join~{1./15};

y1=y$arecp*(y$g1p-y$aH*y$g1)*y$jl;
y2=y$arecp*(y$g2p-y$aH*y$g2)*y$jl;
val1=((y1[[1;;-2]]+y1[[2;;-1]])*intervalLengths/2)//Reverse//Accumulate//Reverse;
val2=((y2[[1;;-2]]+y2[[2;;-1]])*intervalLengths/2)//Reverse//Accumulate//Reverse;

fT$tensor=(y$g2[[1;;-2]]*val1-y$g1[[1;;-2]]*val2)~Join~{0};

lmode=Table[i,{i,1,9}]~Join~Table[10i,{i,1,9}]~Join~Table[100i,{i,1,5}];
ans=Table[0,Length[lmode]];

tab2=Table[((l+2)/(2l+1) SphericalBesselJ[l-1,klist[[uu]] (τ0-τrei)]-(l-1)/(2l+1) SphericalBesselJ[l+1,klist[[uu]] (τ0-τrei)])^2/.l->lmode[[i]],{i,1,Length[lmode]}];

tabvr=tabvr[[;;,1;;totpoint]];
tabvi=tabvi[[;;,1;;totpoint]];
tabvrp=tabvrp[[;;,1;;totpoint]];
tabvip=tabvip[[;;,1;;totpoint]];

Do[
wwmin=FirstPosition[klist,SelectFirst[klist,#>=Abs[klist[[uu]]-klist[[vv]]]&]][[1]];wwmax=FirstPosition[klist,SelectFirst[Reverse[klist],#<=(klist[[uu]]+klist[[vv]])&]][[1]];
Do[
intvr=(2/Mpl^2)*y$arecp*fT$tensor*(klist[[vv]]klist[[ww]](tabvr[[vv]]*tabvr[[ww]]-tabvi[[vv]]*tabvi[[ww]])+(tabvrp[[vv]]*tabvrp[[ww]]-tabvip[[vv]]*tabvip[[ww]]));
intvi=(2/Mpl^2)*y$arecp*fT$tensor*(klist[[vv]]klist[[ww]](tabvr[[vv]]*tabvi[[ww]]+tabvi[[vv]]*tabvr[[ww]])+(tabvrp[[vv]]*tabvip[[ww]]+tabvip[[vv]]*tabvrp[[ww]]));
Vr=Total[(intvr[[1;;-2]]+intvr[[2;;-1]])*intervalLengths/2];
Vi=Total[(intvi[[1;;-2]]+intvi[[2;;-1]])*intervalLengths/2];

ans+=((9/π^3)optdep$rei^2klist[[uu]]dk[[uu]]klist[[vv]]dk[[vv]]klist[[ww]]dk[[ww]]Θ[klist[[uu]],klist[[vv]],klist[[ww]]](Vr^2+Vi^2))*tab2;

,{ww,wwmin,wwmax}];
,{vv,1,Ntot}];

Export[dumpdir<>"uu_"<>$CommandLine[[-1]]<>".dat",MapThread[{#1,#2}&,{lmode,ans}]];
*)
