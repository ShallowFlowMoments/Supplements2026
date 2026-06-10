#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Manufactured Solution*)


(* ::Text:: *)
(*For a given solution, an artificial right-hand side\[LongDash]the "force"\[LongDash]is added to the equations. The numerical solution is then compared with the given solution.*)


(* ::Subsection:: *)
(*Define Simulation Functions*)


LinearMDSM0[params_,cfl_,t_,resolution_]:=DSMSolve[varnamesLDSM0,LDSM0eqns,forceLDSM0,dictLDSM0,initSolLinear[4,mus,omegas,ds,ms],constrLDSM0,params,cfl,t,resolution];
SemilinearMDSM0[params_,cfl_,t_,resolution_]:=DSMSolve[varnamesDSM0,DSM0eqns,forceSLDSM0,dictDSM0,initSolNonlinear[4,mus,omegas,ds,ms],constrSLDSM0,params,cfl,t,resolution];
NonlinearMDSM0[params_,cfl_,t_,resolution_]:=DSMSolve[varnamesDSM0,DSM0eqns,forceDSM0,dictDSM0,initSolNonlinear[4,mus,omegas,ds,ms],constrDSM0,params,cfl,t,resolution];

LinearMDSM1[params_,cfl_,t_,resolution_]:=DSMSolve[varnamesLDSM1,LDSM1eqns,forceLDSM1,dictLDSM1,initSolLinear[7,mus,omegas,ds,ms],constrLDSM1,params,cfl,t,resolution];
SemilinearMDSM1[params_,cfl_,t_,resolution_]:=DSMSolve[varnamesDSM1,DSM1eqns,forceSLDSM1,dictDSM1,initSolNonlinear[7,mus,omegas,ds,ms],constrSLDSM1,params,cfl,t,resolution];
NonlinearMDSM1[params_,cfl_,t_,resolution_]:=DSMSolve[varnamesDSM1,DSM1eqns,forceDSM1,dictDSM1,initSolNonlinear[7,mus,omegas,ds,ms],constrDSM1,params,cfl,t,resolution];

LinearMDSM2[params_,cfl_,t_,resolution_]:=DSMSolve[varnamesLDSM2,LDSM2eqns,forceLDSM2,dictLDSM2,initSolLinear[10,mus,omegas,ds,ms],constrLDSM2,params,cfl,t,resolution];
SemilinearMDSM2[params_,cfl_,t_,resolution_]:=DSMSolve[varnamesDSM2,DSM2eqns,forceSLDSM2,dictDSM2,initSolNonlinear[10,mus,omegas,ds,ms],constrSLDSM2,params,cfl,t,resolution];
NonlinearMDSM2[params_,cfl_,t_,resolution_]:=DSMSolve[varnamesDSM2,DSM2eqns,forceDSM2,dictDSM2,initSolNonlinear[10,mus,omegas,ds,ms],constrDSM2,params,cfl,t,resolution];

LinearMDSM3[params_,cfl_,t_,resolution_]:=DSMSolve[varnamesLDSM3,LDSM3eqns,forceLDSM3,dictLDSM3,initSolLinear[13,mus,omegas,ds,ms],constrLDSM3,params,cfl,t,resolution];
SemilinearMDSM3[params_,cfl_,t_,resolution_]:=DSMSolve[varnamesDSM3,DSM3eqns,forceSLDSM3,dictDSM3,initSolNonlinear[13,mus,omegas,ds,ms],constrSLDSM3,params,cfl,t,resolution];
NonlinearMDSM3[params_,cfl_,t_,resolution_]:=DSMSolve[varnamesDSM3,DSM3eqns,forceDSM3,dictDSM3,initSolNonlinear[13,mus,omegas,ds,ms],constrDSM3,params,cfl,t,resolution];


(* ::Subsection:: *)
(*Set Global Parameters*)


bcType="Periodic" (*"Periodic" or "Fixed"*);
x1=-1.0;
x2=1.0;
xOffset=0;

bshape[x_]= \[Eta] Exp[-x^2/(2 \[Sigma]^2)];
Hb[x_]=bshape[x]-bshape[1];

(*omegas= 3/2{1,1,1,1}\[Pi];
ds={1,1,0,0};
ms={2,1(*1*),2,2(*1*) };
mus={0.1,.1,-0.1,+0.1};
*)

omegas= 2/3{1,1,1,1}\[Pi];
ds={1,1,0,0};
ms={2,1(*1*),2,1(*1*) };
mus={0.1,.1,-0.1,0.1};


(* ::Subsection:: *)
(*Load/define Helper Functions *)


Get[NotebookDirectory[]<>#]&/@{"constraints.wl",
"equations.wl","dictionaries.wl","DSMSolve.wl","visual.wl","sol.wl","forces.wl","varnames.wl"};


ToPureVars[ncUraw_]:=Join[{ncUraw[[1]]},1/ncUraw[[1]]#&/@ncUraw[[2;;-Quotient[Length@ncUraw,3]-1]],ncUraw[[-Quotient[Length@ncUraw,3];;]]]
ToDimlessVars[ncUraw_]:=If[2==Length@ncUraw (*nlsw*),ncUraw/{h0,h0 u0},
(*else: dsm or sme*)
Join[
{ncUraw[[1]]/h0},
1/(h0 u0)#&/@ncUraw[[2;;2+Quotient[Length@ncUraw,3]-1]],
1/(h0 h0 u0)#&/@ncUraw[[2+Quotient[Length@ncUraw,3];;2+2Quotient[Length@ncUraw,3]-1]],
1/(h0)ncUraw[[-Quotient[Length@ncUraw,3];;]]]];
ToDimlessVarsLinear[ncUraw_]:=If[2==Length@ncUraw (*nlsw*),ncUraw/{h0, u0},
(*else: dsm or sme*)
Join[
{ncUraw[[1]]/h0},
1/( u0)#&/@ncUraw[[2;;2+Quotient[Length@ncUraw,3]-1]],
1/(h0  u0)#&/@ncUraw[[2+Quotient[Length@ncUraw,3];;2+2Quotient[Length@ncUraw,3]-1]],
1/(h0)ncUraw[[-Quotient[Length@ncUraw,3];;]]]];

L1ErrH[]:=Total[Abs[Uraw[[1]]-hinit[pts,time]]]/nx;


(* ::Subsection:: *)
(*Visual Output*)


show=True; 
Button[" Show/Hide ",show=!show,BaseStyle->{"GenericButton",12}]
Dynamic[
If[show,
If[Length[Uraw]>0,
(*Visual[ToPureVars@ToDimlessVars@Uraw,400]*)
Visual[Uraw,400],
"no data"
],""],(*,SynchronousUpdating\[Rule]True*)]


h0=0.08;u0=.1;Fr=u0/Sqrt[h0];H=h0/2;Tend= 0.4 0.05;tend= Tend  2/ u0;\[Eta]=0; (*no bathymetry for manufactured sol*)
\[Sigma]=.15;


(* ::Subsection:: *)
(*Manufactured Solutions Linear*)


UrawL0Force=DSMSolve[varnamesLDSM0,LDSM0eqns,forceLDSM0,dictLDSM0,initSolLinear[4,mus,omegas,ds,ms],constrLDSM0,loadParameters,.9, .5,80];


UrawL1Force=DSMSolve[varnamesLDSM1,LDSM1eqns,forceLDSM1,dictLDSM1,initSolLinear[7,mus,omegas,ds,ms],constrLDSM1,loadParameters,.9,Tend,180];


UrawL2Force=DSMSolve[varnamesLDSM2,LDSM2eqns,forceLDSM2,dictLDSM2,initSolLinear[10,mus,omegas,ds,ms],constrLDSM2,loadParameters,.9,Tend,80];


UrawL3Force=DSMSolve[varnamesLDSM3,LDSM3eqns,forceLDSM3,dictLDSM3,initSolLinear[13,mus,omegas,ds,ms],constrLDSM3,loadParameters,.9,Tend,80];


(* ::Subsection:: *)
(*Manufactured Solutions Nonlinear*)


Uraw0Force=DSMSolve[varnamesDSM0,DSM0eqns,forceDSM0,dictDSM0,initSolNonlinear[4,mus,omegas,ds,ms],constrDSM0,loadParameters,.9,Tend,180];


Uraw1Force=DSMSolve[varnamesDSM1,DSM1eqns,forceDSM1,dictDSM1,initSolNonlinear[7,mus,omegas,ds,ms],constrDSM1,loadParameters,.9,Tend,80];


Uraw2Force=DSMSolve[varnamesDSM2,DSM2eqns,forceDSM2,dictDSM2,initSolNonlinear[10,mus,omegas,ds,ms],constrDSM2,loadParameters,.9,Tend,80];


Uraw3Force=DSMSolve[varnamesDSM3,DSM3eqns,forceDSM3,dictDSM3,initSolNonlinear[13,mus,omegas,ds,ms],constrDSM3,loadParameters,.9,Tend,80];


(* ::Subsection:: *)
(*Linear: 2nd order convergence*)


m=80;nsim=3(*4*);


nxVals=Table[m i, {i,1,nsim}]
ref1=Table[{m i,0.01 1./(m i)}, {i,1,nsim}] ;
ref2=Table[{m i,.25 1./(m i)^2}, {i,1,nsim}] ;
ref3=Table[{m i,2.75 1./(m i)^3}, {i,1,nsim}] ;
errs0={#,(LinearMDSM0[loadParameters,.9,tend,#];L1ErrH[])}&/@nxVals


nxVals=Table[m i, {i,1,nsim}]
ref1=Table[{m i,0.01 1./(m i)}, {i,1,nsim}] ;
ref2=Table[{m i,.25 1./(m i)^2}, {i,1,nsim}] ;
ref3=Table[{m i,2.75 1./(m i)^3}, {i,1,nsim}] ;
errs1={#,(LinearMDSM1[loadParameters,.9,tend,#];L1ErrH[])}&/@nxVals


nxVals=Table[m i, {i,1,nsim}]
ref1=Table[{m i,0.01 1./(m i)}, {i,1,nsim}] ;
ref2=Table[{m i,.25 1./(m i)^2}, {i,1,nsim}] ;
ref3=Table[{m i,2.75 1./(m i)^3}, {i,1,nsim}] ;
errs2={#,(LinearMDSM2[loadParameters,.9,tend,#];L1ErrH[])}&/@nxVals


nxVals=Table[m i, {i,1,nsim}]
ref1=Table[{m i,0.01 1./(m i)}, {i,1,nsim}] ;
ref2=Table[{m i,.25 1./(m i)^2}, {i,1,nsim}] ;
ref3=Table[{m i,2.75 1./(m i)^3}, {i,1,nsim}] ;
errs3={#,(LinearMDSM3[loadParameters,.9,tend,#];L1ErrH[])}&/@nxVals


ref1=Table[{m i,   0.003 /(m i)}, {i,1,nsim}] ;
ref2=Table[{m i, .08 1./(m i)^2}, {i,1,nsim}] ;
ref3=Table[{m i,1.75 1./(m i)^3}, {i,1,nsim}] ;


refPlot=ListLogLogPlot[{ref1,ref2},Joined->True,PlotLegends->{"Reference O(nx^(-1))", "Reference O(nx^(-2)"},PlotStyle->{ColorData[97,3],ColorData[97,4]}];
dataPlot=ListLogLogPlot[{errs0,errs1,errs2,errs3},ImageSize->Medium,PlotMarkers->{Automatic, Medium},PlotLegends->{"Data"}(*,ImagePadding->{{50,20},{30,20}}*),Frame->True,FrameStyle->
Directive[Black,Thick],ImageSize->Medium,
LabelStyle->Directive[Black, FontSize->22],PlotRange->All,
GridLines->Automatic,ImageSize->Medium];Show[dataPlot,refPlot]


c=ColorData[97, "ColorList"][[1;;4]];
L=Framed[PointLegend[c,{"LDSM-0","-1","-2","-3"},LegendMarkers->Automatic],Background->White];
lplot=Show[{dataPlot,refPlot},Epilog->Inset[L,Scaled[{.2,.3}]]]


(* ::Subsection:: *)
(*Nonlinear: 2nd order convergence *)


m=80;nsim=3;
nxVals=Table[m i, {i,1,nsim}]



nlerrs0={#,(NonlinearMDSM0[loadParameters,.9,tend,#];L1ErrH[])}&/@nxVals


(*errs0={{50,0.00015463819863096472`},{100,0.000016611591818311602`},{150,6.148991337310515`*^-6},{200,3.3281176784200757`*^-6}}*)


nxVals=Table[m i, {i,1,nsim}]
ref1=Table[{m i,0.01 1./(m i)}, {i,1,nsim}] ;
ref2=Table[{m i,.25 1./(m i)^2}, {i,1,nsim}] ;
ref3=Table[{m i,2.75 1./(m i)^3}, {i,1,nsim}] ;
nlerrs1={#,(NonlinearMDSM1[loadParameters,.9,tend,#];L1ErrH[])}&/@nxVals


(*errs1={{50,0.00022084970208538168`},{100,0.00001612064042020543`},{150,5.384242648618643`*^-6},{200,2.810383714858828`*^-6}}*)


nxVals=Table[m i, {i,1,nsim}]
ref1=Table[{m i,0.01 1./(m i)}, {i,1,nsim}] ;
ref2=Table[{m i,.25 1./(m i)^2}, {i,1,nsim}] ;
ref3=Table[{m i,2.75 1./(m i)^3}, {i,1,nsim}] ;
nlerrs2={#,(NonlinearMDSM2[loadParameters,.9,tend,#];L1ErrH[])}&/@nxVals


(*errs2={{50,0.00023764599505524592`},{100,0.000015646866159264306`},{150,5.194866252932635`*^-6},{200,2.6069352326390206`*^-6}}*)


nxVals=Table[m i, {i,1,nsim}]
ref1=Table[{m i,0.01 1./(m i)}, {i,1,nsim}] ;
ref2=Table[{m i,.25 1./(m i)^2}, {i,1,nsim}] ;
ref3=Table[{m i,2.75 1./(m i)^3}, {i,1,nsim}] ;
nlerrs3={#,(NonlinearMDSM3[loadParameters,.9,tend,#];L1ErrH[])}&/@nxVals


(*errs3={{50,0.00024201880214802762`},{100,0.000015511580901715027`},{150,5.160819850709338`*^-6},{200,2.512022262607487`*^-6}}*)


ref1=Table[{m i,0.015 /(m i)}, {i,1,nsim}] ;
ref2=Table[{m i, 1.35/(m i)^2}, {i,1,nsim}] ;
ref3=Table[{m i,1.2 /(m i)^3}, {i,1,nsim}] ;


refPlot=ListLogLogPlot[{ref1,ref2},Joined->True,PlotLegends->{"Reference O(nx^(-1))", "Reference O(nx^(-2)"}];
dataPlot=ListLogLogPlot[{nlerrs0,nlerrs1,nlerrs2,nlerrs3},ImageSize->Medium,PlotMarkers->{Automatic, Medium},PlotLegends->{"Data"},
(*ImagePadding->{{50,20},{30,20}},*)Frame->True,FrameStyle->
Directive[Black,Thick],ImageSize->Medium,
LabelStyle->Directive[Black, FontSize->22],PlotStyle->{ColorData[97,3],ColorData[97,4]},PlotRange->All,
GridLines->Automatic,ImageSize->Medium];Show[dataPlot,refPlot];


c=ColorData[97, "ColorList"][[1;;4]];L=Framed[PointLegend[c,{"DSM-0","-1","-2","-3"},LegendMarkers->Automatic],Background->White];
Show[{dataPlot,refPlot},ImageSize->Medium,Epilog->Inset[L,Scaled[{.2,.3}]]]


(* export data *)
(*
Export[NotebookDirectory[]<>"manufactured/errs0.csv",errs0]
Export[NotebookDirectory[]<>"manufactured/errs1.csv",errs1]
Export[NotebookDirectory[]<>"manufactured/errs2.csv",errs2]
Export[NotebookDirectory[]<>"manufactured/errs3.csv",errs3]

Export[NotebookDirectory[]<>"manufactured/nlerrs0.csv",nlerrs0]
Export[NotebookDirectory[]<>"manufactured/nlerrs1.csv",nlerrs1]
Export[NotebookDirectory[]<>"manufactured/nlerrs2.csv",nlerrs2]
Export[NotebookDirectory[]<>"manufactured/nlerrs3.csv",nlerrs3]
*)
