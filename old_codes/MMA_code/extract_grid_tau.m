workdir="/tomerv/tomerv-shared/Sida/b-mode-log-griding/f_"<>$CommandLine[[-5]]<>"_L_"<>$CommandLine[[-4]]<>"_tau_"<>$CommandLine[[-3]]<>"_N_"<>$CommandLine[[-2]]<>"_kmax_"<>$CommandLine[[-1]]<>"_nat_tau/parameters/";

If[DirectoryQ[workdir]==False,CreateDirectory[workdir]];

Export[workdir<>"xgrid_HD.dat",#[[1]]&/@Import[workdir<>"../ndsolve/phi.dat","Table"]];
Export[workdir<>"xgrid_cpp.dat",(10^(-40)*#[[1]])&/@Import[workdir<>"../ndsolve/phi.dat","Table"]];

{klist,dk,α,Λ,f,m,τosc}=Import[workdir<>"../parameters.mx"];

klist=1.0*klist[[2;;-1]];
dk=1.0*dk[[2;;-1]];

Export[workdir<>"klist.dat",klist];
Export[workdir<>"dk.dat",dk];

Export[workdir<>"klist_cpp.dat",10^(40)*klist];
Export[workdir<>"dk_cpp.dat",10^(40)*dk];

klist=SetPrecision[klist,100];
dk=SetPrecision[dk,100];

Export[workdir<>"params.dat",{α,Λ,f,m,τosc}];

(*Export[workdir<>"xgrid_trunc.dat",#[[1]]&/@Import[workdir<>"ndsolve_trunc/phi.dat","Table"]];*)
