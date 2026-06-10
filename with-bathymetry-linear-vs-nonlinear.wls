#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Flow Over Uneven Bottom Topography: Linear vs. Nonlinear Equations*)


(* ::Text:: *)
(*This experiment takes place in a parameter set that can only roughly be seen as a perturbation of a stationary flow. Therefore the linear and the nonlinear simulation results diverge.*)


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


\[Sigma]=20/100;\[Eta]= 0.3;
Tend=0.25;g=1;Fr=.75;H=18/100;

h0=2 H; u0=Fr Sqrt[g h0];tend=1/2 Tend 2/ u0 


(* ::Subsection:: *)
(*Compare Linear Systems to the Nonlinear Equivalent*)


(* ::Subsubsection:: *)
(*DSM0*)


(* ::Input:: *)
(*dsm0Nonlinear=ToDimlessVars@DSMSolve[varnamesDSM0,DSM0eqns,ConstantArray[0,Length@varnamesDSM0]&,dictDSM0,initSMEcompare[4,h0,u0,0],constrDSM0,loadParameters,.995,  tend,n];*)
(*dsm0Linear=ToDimlessVarsLinear@DSMSolve[varnamesLDSM0,LDSM0eqns,ConstantArray[0,Length@varnamesLDSM0]&,dictLDSM0,initSMEcompareLinear[4,h0,u0,0],constrLDSM0,loadParameters,.995,tend,n];*)
(*sme0Linear=ToDimlessVarsLinear@DSMSolve[varnamesLDSM0,LDSM0eqns,ConstantArray[0,Length@varnamesLDSM0]&,dictLDSM0,initSMEcompareLinear[4,h0,u0,0],{},loadParameters,.995,tend,n];*)


Grid[{{
(*h*)
Window[{{(dsm0Linear)[[1]]},{(dsm0Nonlinear)[[1]]},{(sme0Linear)[[1]]}},"dsm0: h",400,{0.6,1.2}]},{
(*um*)
Window[{{(dsm0Linear)[[2]]},{(dsm0Nonlinear)[[2]]/dsm0Nonlinear[[1]]},{(sme0Linear)[[2]]}+1},"dsm0: um",400,{.7,1.3}]
(*alpha1*),
Window[{},"dsm0: alpha1",400,{-1,1}]},{
(*wm*)
Window[{{(dsm0Linear)[[3]]},{(dsm0Nonlinear)[[3]]},{(sme0Linear)[[3]]}},"dsm0: wm",400,{-2,2}]
(*gamma1*),
Window[{},"dsm0: gamma1",400,{-.8,.8}]},{
(*q*)
Window[{{(dsm0Linear)[[4]]},{(dsm0Nonlinear)[[4]]},{(sme0Linear)[[4]]}},"dsm0: qm",400,{-.3,.3}]
(*kappa1*),
Window[{},"dsm0: kappa1",400,{-.05,.05}]

}}]


(* ::Subsubsection:: *)
(*DSM1*)


(* ::Input:: *)
(*dsm1Linear=ToDimlessVarsLinear@DSMSolve[varnamesLDSM1,LDSM1eqns,ConstantArray[0,Length@varnamesDSM1]&,dictLDSM1,initSMEcompareLinear[7,h0,u0,0],constrLDSM1,loadParameters,.995,  tend,n];*)
(*sme1Linear=ToDimlessVarsLinear@DSMSolve[varnamesLDSM1,LDSM1eqns,ConstantArray[0,Length@varnamesDSM1]&,dictLDSM1,initSMEcompareLinear[7,h0,u0,0],{},loadParameters,.995,  tend,n];*)
(*dsm1Nonlinear=ToDimlessVars@DSMSolve[varnamesDSM1,DSM1eqns,ConstantArray[0,Length@varnamesDSM1]&,dictDSM1,initSMEcompare[7,h0,u0,0],constrDSM1,loadParameters,.995,  tend,n];*)


sme1Linear //TableForm


(* ::Input:: *)
(*Grid[{{*)
(*(*h*)*)
(*Window[{{(dsm1Linear)[[1]]},{(ToPureVars@dsm1Nonlinear)[[1]]},{(sme1Linear)[[1]]}},"dsm1: h",400,{0.6,1.2}]},{*)
(*(*um*)*)
(*Window[{{(dsm1Linear)[[2]]},{(ToPureVars@dsm1Nonlinear)[[2]]},{(sme1Linear)[[2]]}},"dsm1: um",400,{0.7,1.3}]*)
(*(*alpha1*),*)
(*Window[{{(dsm1Linear)[[3]]},{(ToPureVars@dsm1Nonlinear)[[3]]},{(sme1Linear)[[3]]}},"dsm1: alpha1",400,{-0.5,.5}]},{*)
(*Window[{{(dsm1Linear)[[4]]},{(ToPureVars@dsm1Nonlinear)[[4]]},{(sme1Linear)[[4]]}},"dsm1: wm",400,{-2,2}]*)
(*(*gamma1*),*)
(*Window[{{(dsm1Linear)[[5]]},{(ToPureVars@dsm1Nonlinear)[[5]]},{(sme1Linear)[[5]]}},"dsm1: gamma1",400,{-2,2}]},{*)
(*(*q*)*)
(*Window[{{(dsm1Linear)[[6]]},{(ToPureVars@dsm1Nonlinear)[[6]]},{(sme1Linear)[[6]]}},"dsm1: qm",400,{-.5,.5}]*)
(*(*kappa1*),*)
(*Window[{{(dsm1Linear)[[7]]},{(ToPureVars@dsm1Nonlinear)[[7]]},{(sme1Linear)[[7]]}},"dsm1: kappa1",400,{-.5,.5}]*)
(**)
(*}}]*)
(**)


(*Export[NotebookDirectory[]<>"linear-vs-nonlinear/lsme1.csv",sme1Linear]
Export[NotebookDirectory[]<>"linear-vs-nonlinear/ldsm1.csv",dsm1Linear]
Export[NotebookDirectory[]<>"linear-vs-nonlinear/dsm1.csv",ToPureVars@dsm1Nonlinear]*)
