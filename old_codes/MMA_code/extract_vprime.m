workdir="/tomerv3/Sida/b-mode-log-griding/f_"<>$CommandLine[[-6]]<>"_L_"<>$CommandLine[[-5]]<>"_tau_"<>$CommandLine[[-4]]<>"_N_"<>$CommandLine[[-3]]<>"_kmax_"<>$CommandLine[[-2]]<>"_nat_tau/ndsolve/";

i=ToExpression[$CommandLine[[-1]]];
xgrid=Import[workdir<>"../parameters/xgrid_HD.dat","Table"]//Flatten;

vr=Interpolation[Import[workdir<>"vrpp"<>ToString[i]<>".dat","Table"]];
vi=Interpolation[Import[workdir<>"vipp"<>ToString[i]<>".dat","Table"]];
Export[workdir<>"vrppp"<>ToString[i]<>".dat",MapThread[{#1,vr'[#1]}&,{xgrid}]];
Export[workdir<>"vippp"<>ToString[i]<>".dat",MapThread[{#1,vi'[#1]}&,{xgrid}]];

