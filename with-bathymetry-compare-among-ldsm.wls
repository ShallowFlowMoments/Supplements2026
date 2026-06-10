#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Flow Over Uneven Bottom Topography: Comparing LDSM systems*)


(* ::Text:: *)
(*Waves caused by an uneven bottom topography are simulated in this experiment. Using the default parameters, the flow deviates only slightly from a stationary flow, which justifies the usage of the linear equations.*)


(* ::Subsection:: *)
(*Set Global Parameters*)


n=140;
bcType="Periodic" (*"Periodic" or "Fixed"*);
x1=-1;
x2=1;
xOffset=0; 
bshape[x_]= \[Eta] Exp[-x^2/(2 \[Sigma]^2)];
Hb[x_]=h0 (bshape[x]-bshape[1]);


(* ::Subsection:: *)
(*Load or Define Helper Functions *)


Get[NotebookDirectory[]<>#]&/@{"constraints.wl",
"equations.wl","dictionaries.wl","DSMSolve.wl","visual.wl","sol.wl","forces.wl","varnames.wl"};

ToPureVars[ncUraw_]:=Join[{ncUraw[[1]]},1/ncUraw[[1]]#&/@ncUraw[[2;;-Quotient[Length@ncUraw,3]-1]],ncUraw[[-Quotient[Length@ncUraw,3];;]]];
ToDimlessVars[ncUraw_]:=If[2==Length@ncUraw (*nlsw*),ncUraw/{h0,h0 u0},
(*else: dsm or sme*)
Join[
{ncUraw[[1]]/h0},
1/(h0 u0)#&/@ncUraw[[2;;2+Quotient[Length@ncUraw,3]-1]],
1/(h0 h0/2 u0)#&/@ncUraw[[2+Quotient[Length@ncUraw,3];;2+2Quotient[Length@ncUraw,3]-1]],
1/(h0)ncUraw[[-Quotient[Length@ncUraw,3];;]]]];

ToDimlessVarsLinear[ncUraw_]:=If[2==Length@ncUraw (*nlsw*),ncUraw/{h0, u0},
(*else: dsm or sme*)
Join[
{ncUraw[[1]]/h0},
1/( u0)#&/@ncUraw[[2;;2+Quotient[Length@ncUraw,3]-1]],
1/( h0/2 u0)#&/@ncUraw[[2+Quotient[Length@ncUraw,3];;2+2Quotient[Length@ncUraw,3]-1]],
1/h0 ncUraw[[-Quotient[Length@ncUraw,3];;]]]]+Join[{1,1},ConstantArray[0,Length@ncUraw-2]];

Windowh[UList_,dispname_,iSize_,range_]:=
Show[Plot[{},{x,x1,x2},LabelStyle->Directive[Black,FontSize->22],
FrameStyle->Directive[Black,Thick],ImageSize->Medium,PlotLabel->If[""!=dispname,dispname<>"(x,t)",None],
PlotRange->{{x1,x2},range},ImagePadding->{{100,20},{30,If[""!=dispname,20,10]}},Frame->True,ImageSize->iSize],
MapThread[ListPlot[Transpose[{pts,#1[[1]]}],Joined->True,
PlotStyle->ColorData[106,#2],ImageSize->iSize,PlotRange->All]&,{UList,Range[Length[UList]]}],Background->White,GridLines->Automatic,GridLinesStyle->Black]
Windowhbh[UList_,dispname_,iSize_,range_]:=Show[
Plot[{},{x,x1,x2},LabelStyle->Directive[Black,FontSize->22],FrameStyle->Directive[Black,Thick]
,ImageSize->Medium,PlotLabel->If[""!=dispname,dispname<>"(x,t)",None],PlotRange->{{x1,x2},range},ImagePadding->{{100,20},{30,If[""!=dispname,20,10]}},Frame->True,ImageSize->iSize],
MapThread[ListPlot[Transpose[{pts,#1[[1]]}],Joined->True,Filling->0,
FillingStyle->Directive[RGBColor[0.4,0.6,1],Opacity[.3]],PlotStyle->ColorData[106,#2],ImageSize->iSize,PlotRange->All]&,{UList,Range[Length[UList]]}],Background->White,GridLines->None,GridLinesStyle->Black]
Windowhb[UList_,dispname_,iSize_,range_]:=Show[
Plot[{},{x,x1,x2},LabelStyle->Directive[Black,FontSize->22],FrameStyle->Directive[Black,Thick]
,ImageSize->Medium,PlotLabel->If[""!=dispname,dispname<>"(x,t)",None],PlotRange->{{x1,x2},range},
ImagePadding->{{100,20},{30,If[""!=dispname,20,10]}},Frame->True,ImageSize->iSize],
MapThread[ListPlot[Transpose[{pts,#1[[1]]}],Joined->True,PlotStyle->Darker@LightBrown,Filling->-.025,FillingStyle->Opacity[1],ImageSize->iSize,
PlotRange->All]&,{UList,Range[Length[UList]]}],Background->White,GridLines->None,GridLinesStyle->Black]


(* ::Subsection:: *)
(*Visual Output*)


(* ::Input::Initialization:: *)
show=True; 
Button[" Show/Hide ",show=!show,BaseStyle->{"GenericButton",12}]
Dynamic[
If[show,
If[Length[Uraw]>0,
VisualList[{{ToDimlessVarsLinear[Uraw][[1]]+HbList/h0,ToDimlessVarsLinear[Uraw][[2]]}},{"h","hu"},400,{{.9,1.1},{.9,1.1}}],
"no data"
],""],(*,SynchronousUpdating\[Rule]True*)]


(* ::Subsection:: *)
(*Set Experiment Parameters*)


\[Sigma]=20/100;\[Eta]= 0.03;
Tend=0.25;g=1;Fr=.75;H=8/100;

h0=2 H; u0=Fr Sqrt[g h0];tend=1/2 Tend 2/ u0;


(* ::Subsection:: *)
(*Compare Linear Systems by Level*)


sme0Linear=ToDimlessVarsLinear@DSMSolve[varnamesLDSM0,LDSM0eqns,ConstantArray[0,Length@varnamesLDSM0]&,dictLDSM0,initSMEcompareLinear[4,h0,u0,0],{},loadParameters,.995,tend,n];


sme3Linear=ToDimlessVarsLinear@DSMSolve[varnamesDSM3,LDSM3eqns,ConstantArray[0,Length@varnamesDSM3]&,dictLDSM3,initSMEcompareLinear[13,h0,u0,0],{},loadParameters,.995,  tend,n];


dsm0Linear=ToDimlessVarsLinear@DSMSolve[varnamesLDSM0,LDSM0eqns,ConstantArray[0,Length@varnamesLDSM0]&,dictLDSM0,initSMEcompareLinear[4,h0,u0,0],constrLDSM0,loadParameters,.995,tend,n];


dsm1Linear=ToDimlessVarsLinear@DSMSolve[varnamesLDSM1,LDSM1eqns,ConstantArray[0,Length@varnamesDSM1]&,dictLDSM1,initSMEcompareLinear[7,h0,u0,0],constrLDSM1,loadParameters,.995,  tend,n];


dsm2Linear=ToDimlessVarsLinear@DSMSolve[varnamesDSM2,LDSM2eqns,ConstantArray[0,Length@varnamesDSM2]&,dictLDSM2,initSMEcompareLinear[10,h0,u0,0],constrLDSM2,loadParameters,.995,  tend,n];


dsm3Linear=ToDimlessVarsLinear@DSMSolve[varnamesDSM3,LDSM3eqns,ConstantArray[0,Length@varnamesDSM3]&,dictLDSM3,initSMEcompareLinear[13,h0,u0,0],constrLDSM3,loadParameters,.995,  tend,n];


hbhplot=Show[{Windowhbh[{{(dsm0Linear)[[1]]+HbList/h0},{(dsm1Linear)[[1]]+HbList/h0},{(dsm2Linear)[[1]]+HbList/h0},
{(dsm3Linear)[[1]]+HbList/h0},{(sme3Linear)[[1]]+HbList/h0}},"",400,{-.025,1.2}]
,
Windowhb[{{HbList/h0}},"",400,{0,1.2}]
}]


hplot=Windowh[{{(dsm0Linear)[[1]]},{(dsm1Linear)[[1]]},{(dsm2Linear)[[1]]},{(dsm3Linear)[[1]]},{(sme3Linear)[[1]]}},"",400,{0.95,1.02}]
(*um*)
umplot=Windowh[{{(dsm0Linear)[[2]]},{(dsm1Linear)[[2]]},{(dsm2Linear)[[2]]},{(dsm3Linear)[[2]]},{(sme3Linear)[[2]]}},"",400,{.95,1.05}]
(*alpha1*)
a1plot=Windowh[{{(dsm1Linear)[[3]]},{(dsm1Linear)[[3]]},{(dsm2Linear)[[3]]},{(dsm3Linear)[[3]]},{sme3Linear[[3]]}},"",400,{-0.05,.05}]
(*wm*)
wmplot=Windowh[{{(dsm0Linear)[[3]]},{(dsm1Linear)[[4]]},{(dsm2Linear)[[5]]},{(dsm3Linear)[[6]]}},"",400,{-1,1}]
(*gamma1*)
g1plot=Windowh[{{(dsm1Linear)[[5]]},{(dsm1Linear)[[5]]},{(dsm2Linear)[[6]]},{(dsm3Linear)[[7]]}},"",400,{-1,1}]
(*q*)
qmplot=Windowh[{{(dsm0Linear)[[4]]},{(dsm1Linear)[[6]]},{(dsm2Linear)[[8]]},{(dsm3Linear)[[10]]}},"",400,{-.02,.02}]
(*kappa1*)
k1plot=Windowh[{{(dsm1Linear)[[7]]},{(dsm1Linear)[[7]]},{(dsm2Linear)[[9]]},{(dsm3Linear)[[11]]}},"",400,{-.02,.02}]



(* export data *)


hvals=Flatten[{{(dsm0Linear)[[1]]},{(dsm1Linear)[[1]]},{(dsm2Linear)[[1]]},{(dsm3Linear)[[1]]},{(sme3Linear)[[1]]}},1];


umvals=Flatten[{{(dsm0Linear)[[1]]},{(dsm1Linear)[[2]]},{(dsm2Linear)[[2]]},{(dsm3Linear)[[2]]},{(sme3Linear)[[2]]}},1];


a1vals=Flatten[{{(dsm1Linear)[[3]]},{(dsm2Linear)[[3]]},{(dsm3Linear)[[3]]},{sme3Linear[[3]]}},1];


wmvals=Flatten[{{(dsm0Linear)[[3]]},{(dsm1Linear)[[4]]},{(dsm2Linear)[[5]]},{(dsm3Linear)[[6]]}},1];


g1vals=Flatten[{{(dsm1Linear)[[5]]},{(dsm2Linear)[[6]]},{(dsm3Linear)[[7]]}},1];


qmvals=Flatten[{{(dsm0Linear)[[4]]},{(dsm1Linear)[[6]]},{(dsm2Linear)[[8]]},{(dsm3Linear)[[10]]}},1];


k1vals=Flatten[{{(dsm1Linear)[[7]]},{(dsm2Linear)[[9]]},{(dsm3Linear)[[11]]}},1];


(* Export[NotebookDirectory[]<>"lin-with-bathymetry/h.csv",hvals]
Export[NotebookDirectory[]<>"lin-with-bathymetry/um.csv",umvals]
Export[NotebookDirectory[]<>"lin-with-bathymetry/a1.csv",a1vals]
Export[NotebookDirectory[]<>"lin-with-bathymetry/wm.csv",wmvals]
Export[NotebookDirectory[]<>"lin-with-bathymetry/g1.csv",g1vals]
Export[NotebookDirectory[]<>"lin-with-bathymetry/qm.csv",qmvals]
Export[NotebookDirectory[]<>"lin-with-bathymetry/k1.csv",k1vals] *)


(* ::Subsection:: *)
(*Compare Some Linear Systems to the Nonlinear Equivalent*)


(* ::Subsubsection:: *)
(*DSM0*)


(* ::Input:: *)
(*dsm0Nonlinear=ToDimlessVars@DSMSolve[varnamesDSM0,DSM0eqns,ConstantArray[0,Length@varnamesDSM0]&,dictDSM0,initSMEcompare[4,h0,u0,0],constrDSM0,loadParameters,.995,  tend,n];*)
(*dsm0Linear=ToDimlessVarsLinear@DSMSolve[varnamesLDSM0,LDSM0eqns,ConstantArray[0,Length@varnamesLDSM0]&,dictLDSM0,initSMEcompareLinear[4,h0,u0,0],constrLDSM0,loadParameters,.995,tend,n];*)
(*sme0Linear=ToDimlessVarsLinear@DSMSolve[varnamesLDSM0,LDSM0eqns,ConstantArray[0,Length@varnamesLDSM0]&,dictLDSM0,initSMEcompareLinear[4,h0,u0,0],{},loadParameters,.995,tend,n];*)


Grid[{{
(*h*)
Window[{{(dsm0Linear)[[1]]},{(dsm0Nonlinear)[[1]]},{(sme0Linear)[[1]]}},"dsm0: h",400,{0.95,1.02}]},{
(*um*)
Window[{{(dsm0Linear)[[2]]},{(dsm0Nonlinear)[[2]]/dsm0Nonlinear[[1]]},{(sme0Linear)[[2]]}+1},"dsm0: um",400,{.9,1.2}]
(*alpha1*),
Window[{},"dsm0: alpha1",400,{-0.1,.1}]},{
(*wm*)
Window[{{(dsm0Linear)[[3]]},{(dsm0Nonlinear)[[3]]},{(sme0Linear)[[3]]}},"dsm0: wm",400,{-.7,.7}]
(*gamma1*),
Window[{},"dsm0: gamma1",400,{-.3,.3}]},{
(*q*)
Window[{{(dsm0Linear)[[4]]},{(dsm0Nonlinear)[[4]]},{(sme0Linear)[[4]]}},"dsm0: qm",400,{-.02,.02}]
(*kappa1*),
Window[{},"dsm0: kappa1",400,{-.005,.005}]

}}]


(* ::Subsubsection:: *)
(*DSM1*)


(* ::Input:: *)
(*dsm1Linear=ToDimlessVarsLinear@DSMSolve[varnamesLDSM1,LDSM1eqns,ConstantArray[0,Length@varnamesDSM1]&,dictLDSM1,initSMEcompareLinear[7,h0,u0,0],constrLDSM1,loadParameters,.995,  tend,n];*)
(*sme1Linear=ToDimlessVarsLinear@DSMSolve[varnamesLDSM1,LDSM1eqns,ConstantArray[0,Length@varnamesDSM1]&,dictLDSM1,initSMEcompareLinear[7,h0,u0,0],{},loadParameters,.995,  tend,n];*)
(*dsm1Nonlinear=ToDimlessVars@DSMSolve[varnamesDSM1,DSM1eqns,ConstantArray[0,Length@varnamesDSM1]&,dictDSM1,initSMEcompare[7,h0,u0,0],constrDSM1,loadParameters,.995,  tend,n];*)


(* ::Input:: *)
(*Grid[{{*)
(*(*h*)*)
(*Window[{{(dsm1Linear)[[1]]},{(dsm1Nonlinear)[[1]]},{(sme1Linear)[[1]]}},"dsm1: h",400,{0.9,1.1}]},{*)
(*(*um*)*)
(*Window[{{(dsm1Linear)[[2]]},{(ToPureVars@dsm1Nonlinear)[[2]]},{(sme1Linear)[[2]]+1}},"dsm1: um",400,{0.8,1.2}]*)
(*(*alpha1*),*)
(*Window[{{(dsm1Linear)[[3]]},{(ToPureVars@dsm1Nonlinear)[[3]]},{(sme1Linear)[[3]]}},"dsm1: alpha1",400,{-0.1,.1}]},{*)
(*Window[{{(dsm1Linear)[[4]]},{(ToPureVars@dsm1Nonlinear)[[4]]},{(sme1Linear)[[4]]}},"dsm1: wm",400,{-.7,.7}]*)
(*(*gamma1*),*)
(*Window[{{(dsm1Linear)[[5]]},{(ToPureVars@dsm1Nonlinear)[[5]]},{(sme1Linear)[[5]]}},"dsm1: gamma1",400,{-.3,.3}]},{*)
(*(*q*)*)
(*Window[{{(dsm1Linear)[[6]]},{(ToPureVars@dsm1Nonlinear)[[6]]},{(sme1Linear)[[6]]}},"dsm1: qm",400,{-.02,.02}]*)
(*(*kappa1*),*)
(*Window[{{(dsm1Linear)[[7]]},{(ToPureVars@dsm1Nonlinear)[[7]]},{(sme1Linear)[[7]]}},"dsm1: kappa1",400,{-.005,.005}]*)
(**)
(*}}]*)
(**)
