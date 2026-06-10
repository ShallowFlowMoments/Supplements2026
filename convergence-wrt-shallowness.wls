#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Convergence with Respect to Shallowness*)


(* ::Text:: *)
(*DSM can capture dispersive effects. These effects become more prominent as shallowness decreases . *)
(*This script runs the same experiment multiple times at different degrees of shallowness .*)


(* ::Subsection:: *)
(*Set Global Parameters*)


bcType="Periodic" (*"Periodic" or "Fixed"*);
x1=-1.0;
x2=1.0;
xOffset=0;

n=140;

bshape[x_]= \[Eta] Exp[-x^2/(2 \[Sigma]^2)];
Hb[x_]=bshape[x]-bshape[1];


(* ::Subsection:: *)
(*Load or Define Helper Functions *)


Get[NotebookDirectory[]<>#]&/@{"constraints.wl",
"equations.wl","dictionaries.wl","DSMSolve.wl","visual.wl","sol.wl","forces.wl","varnames.wl"};

ToPureVars[ncUraw_]:=Join[{ncUraw[[1]]},1/ncUraw[[1]]#&/@ncUraw[[2;;-Quotient[Length@ncUraw,3]-1]],ncUraw[[-Quotient[Length@ncUraw,3];;]]];
ToDimlessVarsWaveSpeed[ncUraw_]:=If[2==Length@ncUraw (*nlsw*),ncUraw/{h0,h0 100 Sqrt[h0]},
(*else: dsm or sme*)
Join[
{ncUraw[[1]]/h0},
1/(h0 100 Sqrt[h0])#&/@ncUraw[[2;;2+Quotient[Length@ncUraw,3]-1]],
1/(h0 h0/2 100 Sqrt[h0])#&/@ncUraw[[2+Quotient[Length@ncUraw,3];;2+2Quotient[Length@ncUraw,3]-1]],
1/(h0)ncUraw[[-Quotient[Length@ncUraw,3];;]]]];

Windowh[UList_,dispname_,iSize_,range_]:=
Show[Plot[{},{x,x1,x2},LabelStyle->Directive[Black,FontSize->22],
FrameStyle->Directive[Black,Thick],ImageSize->Medium,PlotLabel->If[""!=dispname,dispname<>"(x,t)",None],
PlotRange->{{x1,x2},range},ImagePadding->{{100,20},{30,If[""!=dispname,20,10]}},Frame->True,ImageSize->iSize],
MapThread[ListPlot[Transpose[{pts,#1[[1]]}],Joined->True,
PlotStyle->ColorData[106,#2],ImageSize->iSize,PlotRange->All]&,{UList,Range[Length[UList]]}],Background->White,GridLines->Automatic,GridLinesStyle->Black]


(* ::Subsection:: *)
(*Visual Output*)


show=True; 
Button[" Show/Hide ",show=!show,BaseStyle->{"GenericButton",12}]
Dynamic[
If[show,
If[Length[Uraw]>0,
VisualList[{ToDimlessVarsWaveSpeed@Uraw},{"h","hu"},400,{{.5,1.5},{.5,1.5}}],
"no data"
],""],(*,SynchronousUpdating\[Rule]True*)]


(* ::Subsection:: *)
(*Set Experiment Parameters*)


g=1;
\[Sigma]=.15;\[Eta]=0.00 (*no bottom topography*);
Fr=.25;Tend=.25;Hvals={0.0025,0.005,0.0075}; (* shallowness parameter has symbol 'H' here *)
s0=0*-u0;tend=Tend 2/u0;


(* ::Subsection:: *)
(*Experiments With Nonlinear Equations*)


sme1=Block[{H=#},
h0=2 H;u0=0 ;tend=2/Sqrt[g h0] (*8.52*);Print@TableForm[{h0,u0,tend}];
ToDimlessVarsWaveSpeed@DSMSolve[varnamesDSM1,DSM1eqns,ConstantArray[0,Length@varnamesDSM1]&,dictDSM1,initSMEcompareNoOffset[7,h0,u0,0],{},loadParameters,.995,tend  ,n]]&/@{Hvals[[1]]};


dsm1pseudoDL=Block[{H=#},
h0=2 H;u0=0 ;tend= 2 / Sqrt[g h0] (*8.52*);Print@TableForm[{h0,u0,tend}];
ToDimlessVarsWaveSpeed@DSMSolve[varnamesDSM1,DSM1eqns,
ConstantArray[0,Length@varnamesDSM1]&,dictDSM1,
initSMEcompareNoOffset[7,h0,u0,0],constrDSM1,loadParameters,.995,  tend,n]]&/@Hvals;


hvals={ToPureVars[#][[1]]}&/@Join[{sme1[[1]]},dsm1pseudoDL[[1;;3]]];
umvals={ToPureVars[#][[2]]}&/@Join[{sme1[[1]]},dsm1pseudoDL[[1;;3]]];
a1vals={ToPureVars[#][[3]]}&/@Join[{sme1[[1]]},dsm1pseudoDL[[1;;3]]];
qmvals=Transpose@{Join[{ConstantArray[0,n]},dsm1pseudoDL[[;;,6]]]};


legend=Framed[SwatchLegend[{RGBColor[0.34398, 0.49112, 0.89936],RGBColor[0.97, 0.606, 0.081],RGBColor[0.91, 0.318, 0.243],RGBColor[0.448, 0.69232, 0.1538]},
{"SME1 (S=0 limit)","DSM1;S=1/400","DSM1;S=1/200  \[ThinSpace]\[VeryThinSpace]","DSM1;S=3/400"},LegendLayout->"Row"],Background->White]


hplot=Windowh[hvals,"h",400,{0.95,1.3}]
umplot=Windowh[umvals,"um",400,{-0.005,0.005}] 
a1plot=Windowh[a1vals,"a1",400,{-0.0002,0.0002}]
qmplot=Windowh[qmvals,"qm",400,{-.02,.02}]


(* export data *)

(* Export[NotebookDirectory[]<>"convergenceH/h.csv",hvals]
Export[NotebookDirectory[]<>"convergenceH/um.csv",umvals]
Export[NotebookDirectory[]<>"convergenceH/a1.csv",a1vals]
Export[NotebookDirectory[]<>"convergenceH/qm.csv",qmvals] *)


(* ::Subsection:: *)
(*Experiments With Linear Equations*)


sme1linear=Block[{H=#},
h0=2 H;u0=0 ;tend=2/Sqrt[g h0] (*8.52*);Print@TableForm[{h0,u0,tend}];
ToDimlessVarsWaveSpeed@DSMSolve[varnamesLDSM1,LDSM1eqns,ConstantArray[0,Length@varnamesDSM1]&,dictLDSM1,initSMEcompareLinearNoOffset[7,h0,u0,0],{},loadParameters,.995,tend  ,n]]&/@{Hvals[[1]]};


ldsm1pseudoDL=Block[{H=#},
h0=2 H;u0=0 ;tend= 2 / Sqrt[g h0] (*8.52*);Print@TableForm[{h0,u0,tend}];
ToDimlessVarsWaveSpeed@DSMSolve[varnamesLDSM1,LDSM1eqns,
ConstantArray[0,Length@varnamesDSM1]&,dictLDSM1,
initSMEcompareLinearNoOffset[7,h0,u0,0],constrLDSM1,loadParameters,.995,  tend,n]]&/@Hvals;


lhvals={Identity[#][[1]]}&/@Join[{sme1linear[[1]]},ldsm1pseudoDL[[1;;3]]];
lumvals={Identity[#][[2]]}&/@Join[{sme1linear[[1]]},ldsm1pseudoDL[[1;;3]]];
la1vals={Identity[#][[3]]}&/@Join[{sme1linear[[1]]},ldsm1pseudoDL[[1;;3]]];
lqmvals=Transpose@{Join[{ConstantArray[0,n]},ldsm1pseudoDL[[;;,6]]]};


legend=Framed[SwatchLegend[{RGBColor[0.34398, 0.49112, 0.89936],RGBColor[0.97, 0.606, 0.081],RGBColor[0.91, 0.318, 0.243],RGBColor[0.448, 0.69232, 0.1538]},
{"LSME1 (S=0 limit)","LDSM1;S=1/400","LDSM1;S=1/200  \[ThinSpace]\[VeryThinSpace]","LDSM1;S=3/400"},LegendLayout->"Row"],Background->White]


lhplot=Windowh[lhvals,"h",400,{0.95,1.3}]
lumplot=Windowh[lumvals,"um",400,{-0.005,0.005}] 
la1plot=Windowh[la1vals,"a1",400,{-0.0002,0.0002}]
lqmplot=Windowh[lqmvals,"q",400,{-.02,.02}]
