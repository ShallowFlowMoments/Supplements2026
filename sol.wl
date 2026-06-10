(* ::Package:: *)

initSolLinear[nVars_,mus_,omegas_,ds_,ms_]:=(
{\[Mu]1,\[Mu]2,\[Mu]3,\[Mu]4}=mus;
{\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4}=omegas;
{d1,d2,d3,d4}=ds;
{m1,m2,m3,m4}=ms;

nMoments=Quotient[nVars,3]-1;

hinit[x_,t_]:=d1+\[Mu]1 Sin[\[Pi] x+\[Omega]1 t^m1](*-Hb[x]*);
	uinit[x_,t_]:=d2+\[Mu]2 Sin[\[Pi] x+\[Omega]2 t^m2];
aMomentsInit[x_,t_]:=Table[\[Mu]2^2 Sin[\[Pi] x+\[Omega]2 t^m2],{i,1,nMoments}];
	winit[x_,t_]:=d3+\[Mu]3 Sin[\[Pi] x+\[Omega]3 t^m3];
gMomentsInit[x_,t_]:=Table[\[Mu]3^2 Sin[\[Pi] x+\[Omega]3 t^m3],{i,1,nMoments}];
	qinit[x_,t_]:=d4+\[Mu]4 Sin[\[Pi] x+ \[Omega]4 t^m4];
kMomentsinit[x_,t_]:=Table[\[Mu]4^2 Sin[\[Pi] x+\[Omega]4 t^m4],{i,1,nMoments}];

Function[{x,t},Flatten[{hinit[x,t],uinit[x,t],aMomentsInit[x,t],winit[x,t],gMomentsInit[x,t],qinit[x,t],kMomentsinit[x,t]}]]
);


initSolNonlinear[nVars_,mus_,omegas_,ds_,ms_]:=(
{\[Mu]1,\[Mu]2,\[Mu]3,\[Mu]4}=mus;
{\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4}=omegas;
{d1,d2,d3,d4}=ds;
{m1,m2,m3,m4}=ms;

nMoments=If[2==nVars,0,Quotient[nVars,3]-1];

hinit[x_,t_]=d1+\[Mu]1 Sin[\[Pi] x+\[Omega]1 t^m1](*-Hb[x]*);
	uinit[x_,t_]=d2+\[Mu]2 Sin[\[Pi] x+\[Omega]2 t^m2];
aMomentsInit[x_,t_]=Table[\[Mu]2^2 Sin[\[Pi] x+\[Omega]2 t^m2],{i,1,nMoments}];
	winit[x_,t_]=d3+\[Mu]3 Sin[\[Pi] x+\[Omega]3 t^m3];
gMomentsInit[x_,t_]=Table[\[Mu]3^2 Sin[\[Pi] x+\[Omega]3 t^m3],{i,1,nMoments}];
	qinit[x_,t_]=d4+\[Mu]4 Sin[\[Pi] x+ \[Omega]4 t^m4];
kMomentsinit[x_,t_]=Table[\[Mu]4^2 Sin[\[Pi] x+\[Omega]4 t^m4],{i,1,nMoments}];

If[2==nVars, Function[{x,t},{hinit[x,t],hinit[x,t]uinit[x,t]}],
Function[{x,t},Flatten[{hinit[x,t],hinit[x,t]uinit[x,t],hinit[x,t]aMomentsInit[x,t],hinit[x,t]winit[x,t],hinit[x,t]gMomentsInit[x,t],qinit[x,t],kMomentsinit[x,t]}]]
]);


initSMEcompareNoOffset[nVars_,h0_,u0_,s0_]:=(

nMoments=If[2==nVars,0,Quotient[nVars,3]-1];
poly[i_,var_]:=LegendreP[i,-2var+1];

hinit[x_,t_]=h0+If[Hb[x]>0,0,1]*h0 Exp[3Cos[\[Pi] (x)]]/Exp[4]-Hb[x];
	uinit[x_,t_]=u0;
aMomentsInit[x_,t_]=Table[If[i==1,s0,0],{i,1,nMoments}];
ustart[x_,\[Zeta]_]=u0+s0 poly[1,\[Zeta]];
hstart[x_]=hinit[x,0];
wstart[x_,\[Zeta]_]= Simplify[-hstart[x] D[Integrate[ustart[xb,\[Zeta]b],{\[Zeta]b,0,\[Zeta]}],xb]+D[hstart[xb],xb](ustart[x,\[Zeta]] \[Zeta]-Integrate[ustart[x,\[Zeta]b],{\[Zeta]b,0,\[Zeta]}]) /.xb->x];
	winit[x_,t_]=If[IntegerQ[wstart[x,\[Zeta]]]&&wstart[x,\[Zeta]]==0,0,Integrate[wstart[x,\[Zeta]],{\[Zeta],0,1}]];
gMomentsInit[x_,t_]=Table[If[IntegerQ[wstart[x,\[Zeta]]]&&wstart[x,\[Zeta]]==0,0,Integrate[wstart[x,\[Zeta]]poly[i,\[Zeta]],{\[Zeta],0,1}]],{i,1,nMoments}];
	qinit[x_,t_]=0;
kMomentsinit[x_,t_]=Table[0,{i,1,nMoments}];

If[2==nVars, Function[{x,t},{hinit[x,t],hinit[x,t]uinit[x,t]}],
Function[{x,t},Flatten[{hinit[x,t],hinit[x,t]uinit[x,t],hinit[x,t]aMomentsInit[x,t],hinit[x,t]winit[x,t],hinit[x,t]gMomentsInit[x,t],qinit[x,t],kMomentsinit[x,t]}]]
]);


initSMEcompareLinearNoOffset[nVars_,h0_,u0_,s0_]:=(

nMoments=If[2==nVars,0,Quotient[nVars,3]-1];
poly[i_,var_]:=LegendreP[i,-2var+1];

hinit[x_,t_]=h0+If[Hb[x]>0,0,1]*h0 Exp[3Cos[\[Pi] (x)]]/Exp[4]-Hb[x];
	uinit[x_,t_]=u0;
aMomentsInit[x_,t_]=Table[If[i==1,s0,0],{i,1,nMoments}];
ustart[x_,\[Zeta]_]=u0+s0 poly[1,\[Zeta]];
hstart[x_]=hinit[x,0];
wstart[x_,\[Zeta]_]= Simplify[-hstart[x] D[Integrate[ustart[xb,\[Zeta]b],{\[Zeta]b,0,\[Zeta]}],xb]+D[hstart[xb],xb](ustart[x,\[Zeta]] \[Zeta]-Integrate[ustart[x,\[Zeta]b],{\[Zeta]b,0,\[Zeta]}]) /.xb->x];
	winit[x_,t_]=If[IntegerQ[wstart[x,\[Zeta]]]&&wstart[x,\[Zeta]]==0,0,Integrate[wstart[x,\[Zeta]],{\[Zeta],0,1}]];
gMomentsInit[x_,t_]=Table[If[IntegerQ[wstart[x,\[Zeta]]]&&wstart[x,\[Zeta]]==0,0,Integrate[wstart[x,\[Zeta]]poly[i,\[Zeta]],{\[Zeta],0,1}]],{i,1,nMoments}];
	qinit[x_,t_]=0;
kMomentsinit[x_,t_]=Table[0,{i,1,nMoments}];

If[2==nVars, Function[{x,t},{hinit[x,t],hinit[x,t]uinit[x,t]}],
Function[{x,t},Flatten[{hinit[x,t],uinit[x,t],aMomentsInit[x,t],winit[x,t],gMomentsInit[x,t],qinit[x,t],kMomentsinit[x,t]}]]
]);


initSMEcompare[nVars_,h0_,u0_,s0_]:=(

nMoments=If[2==nVars,0,Quotient[nVars,3]-1];
poly[i_,var_]:=LegendreP[i,-2var+1];

hinit[x_,t_]=h0+If[Hb[x]>0,0,1]*h0 Exp[3Cos[\[Pi] (x+0.5)]]/Exp[4]-Hb[x];
	uinit[x_,t_]=u0;
aMomentsInit[x_,t_]=Table[If[i==1,s0,0],{i,1,nMoments}];
ustart[x_,\[Zeta]_]=u0+s0 poly[1,\[Zeta]];
hstart[x_]=hinit[x,0];
wstart[x_,\[Zeta]_]= Simplify[-hstart[x] D[Integrate[ustart[xb,\[Zeta]b],{\[Zeta]b,0,\[Zeta]}],xb]+D[hstart[xb],xb](ustart[x,\[Zeta]] \[Zeta]-Integrate[ustart[x,\[Zeta]b],{\[Zeta]b,0,\[Zeta]}]) /.xb->x];
	winit[x_,t_]=(*If[IntegerQ[wstart[x,\[Zeta]]]&&wstart[x,\[Zeta]]==0,0,Integrate[wstart[x,\[Zeta]],{\[Zeta],0,1}]]*)uinit[x,t]Hb'[x]-hinit[x,t]D[uinit[x,t],x];
gMomentsInit[x_,t_]=Table[If[IntegerQ[wstart[x,\[Zeta]]]&&wstart[x,\[Zeta]]==0,0,Integrate[wstart[x,\[Zeta]]poly[i,\[Zeta]],{\[Zeta],0,1}]],{i,1,nMoments}];
	qinit[x_,t_]=0;
kMomentsinit[x_,t_]=Table[0,{i,1,nMoments}];

If[2==nVars, Function[{x,t},{hinit[x,t],hinit[x,t]uinit[x,t]}],
Function[{x,t},Flatten[{hinit[x,t],hinit[x,t]uinit[x,t],hinit[x,t]aMomentsInit[x,t],hinit[x,t]winit[x,t],hinit[x,t]gMomentsInit[x,t],qinit[x,t],kMomentsinit[x,t]}]]
]);


initSMEcompareLinear[nVars_,h0_,u0_,s0_]:=(

nMoments=If[2==nVars,0,Quotient[nVars,3]-1];
poly[i_,var_]:=LegendreP[i,-2var+1];

hinit[x_,t_]=If[Hb[x]>0,0,1]*h0 Exp[3Cos[\[Pi] (x+1/2)]]/Exp[4]-Hb[x];
	uinit[x_,t_]=0;
aMomentsInit[x_,t_]=Table[If[i==1,s0,0],{i,1,nMoments}];
ustart[x_,\[Zeta]_]=s0 poly[1,\[Zeta]];
hstart[x_]=hinit[x,0];
wstart[x_,\[Zeta]_]= Simplify[-hstart[x] D[Integrate[ustart[xb,\[Zeta]b],{\[Zeta]b,0,\[Zeta]}],xb]+D[hstart[xb],xb](ustart[x,\[Zeta]] \[Zeta]-Integrate[ustart[x,\[Zeta]b],{\[Zeta]b,0,\[Zeta]}]) /.xb->x];
	winit[x_,t_]=(*- If[IntegerQ[wstart[x,\[Zeta]]]&&wstart[x,\[Zeta]]==0,0,Integrate[wstart[x,\[Zeta]],{\[Zeta],0,1}]]*)u0 Hb'[x]-h0 D[uinit[t,x],x];
gMomentsInit[x_,t_]=Table[If[IntegerQ[wstart[x,\[Zeta]]]&&wstart[x,\[Zeta]]==0,0,Integrate[wstart[x,\[Zeta]]poly[i,\[Zeta]],{\[Zeta],0,1}]],{i,1,nMoments}];
	qinit[x_,t_]=0;
kMomentsinit[x_,t_]=Table[0,{i,1,nMoments}];


Function[{x,t},Flatten[{hinit[x,t],uinit[x,t],aMomentsInit[x,t],winit[x,t],
gMomentsInit[x,t],qinit[x,t],kMomentsinit[x,t]}]]);




initSolNonlinearDimless[nVars_,mus_,omegas_,ds_,ms_]:=(
{\[Mu]1,\[Mu]2,\[Mu]3,\[Mu]4}=mus;
{\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4}=omegas;
{d1,d2,d3,d4}=ds;
{m1,m2,m3,m4}=ms;

nMoments=If[2==nVars,0,Quotient[nVars,3]-1];

hinit[x_,t_]=1/h0 (d1+\[Mu]1 Sin[\[Pi] x+\[Omega]1 t^m1])(*-Hb[x]*);
	uinit[x_,t_]=1/u0 (d2+\[Mu]2 Sin[\[Pi] x+\[Omega]2 t^m2]);
aMomentsInit[x_,t_]=Table[1/u0 (\[Mu]2^2 Sin[\[Pi] x+\[Omega]2 t^m2]),{i,1,nMoments}];
	winit[x_,t_]=1/(h0 u0)(d3+\[Mu]3 Sin[\[Pi] x+\[Omega]3 t^m3]);
gMomentsInit[x_,t_]=Table[1/(h0 u0)(\[Mu]3^2 Sin[\[Pi] x+\[Omega]3 t^m3]),{i,1,nMoments}];
	qinit[x_,t_]=1/h0 (d4+\[Mu]4 Sin[\[Pi] x+ \[Omega]4 t^m4]);
kMomentsinit[x_,t_]=Table[1/h0 (\[Mu]4^2 Sin[\[Pi] x+\[Omega]4 t^m4]),{i,1,nMoments}];

If[2==nVars, Function[{x,t},{hinit[x,t],hinit[x,t]uinit[x,t]}],
Function[{x,t},Flatten[{hinit[x,t],hinit[x,t]uinit[x,t],hinit[x,t]aMomentsInit[x,t],hinit[x,t]winit[x,t],hinit[x,t]gMomentsInit[x,t],qinit[x,t],kMomentsinit[x,t]}]]
]);
