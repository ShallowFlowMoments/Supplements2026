#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Hierarchical Convergence*)


(* ::Text:: *)
(*The higher the level of a DSM system, the greater the accuracy. This is confirmed by the following code which runs the same experiment multiple times using different DSM levels.*)


(* ::Subsection:: *)
(*Set Global Parameters*)


bcType="Periodic" (*"Periodic" or "Fixed"*);
x1=-1.0;
x2=1.0;
n=140(*280*);
xOffset=0;

bshape[x_]= \[Eta] Exp[-x^2/(2 \[Sigma]^2)];
Hb[x_]=bshape[x]-bshape[1];


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



(* ::Subsection:: *)
(*Visual Output*)


show=True; 
Button[" Show/Hide ",show=!show,BaseStyle->{"GenericButton",12}]
Dynamic[
If[show,
If[Length[Uraw]>0,
VisualList[{ToDimlessVars@Uraw},{"h","hu"},400,{{.5,1.5},{.5,1.5}}],
"no data"
],""],(*,SynchronousUpdating\[Rule]True*)]


(* ::Input:: *)
(**)


(* ::Subsection:: *)
(*Set Experiment Parameters*)


g=1;Fr=.25;H=1/200;Tend=.25;
\[Sigma]=.15;\[Eta]=0.0; (* flat bottom *)
h0=2 H; u0=Fr Sqrt[g h0];tend= Tend 2/ u0; (*dimensional variables for pseudodimensionless simulation*)


(* ::Subsection:: *)
(*Experiments With Nonlinear Equations*)


dsm0[H,n]=DSMSolve[varnamesDSM0,DSM0eqns,ConstantArray[0,Length@varnamesDSM0]&,dictDSM0,initSMEcompare[4,h0,u0,0*-1],constrDSM0,loadParameters,.995,tend,n];
dsm1[H,n]=DSMSolve[varnamesDSM1,DSM1eqns,ConstantArray[0,Length@varnamesDSM1]&,dictDSM1,initSMEcompare[7,h0,u0,0*-1],constrDSM1,loadParameters,.995,tend,n];
dsm2[H,n]=DSMSolve[varnamesDSM2,DSM2eqns,ConstantArray[0,Length@varnamesDSM2]&,dictDSM2,initSMEcompare[10,h0,u0,0*-1],constrDSM2,loadParameters,.995,tend,n];
dsm3[H,n]=DSMSolve[varnamesDSM3,DSM3eqns,ConstantArray[0,Length@varnamesDSM3]&,dictDSM3,initSMEcompare[13,h0,u0,0*-1],constrDSM3,loadParameters,.995,tend,n];   


dsm0dl[H,n]=ToDimlessVars@dsm0[H,n];
dsm1dl[H,n]=ToDimlessVars@dsm1[H,n];
dsm2dl[H,n]=ToDimlessVars@dsm2[H,n];
dsm3dl[H,n]=ToDimlessVars@dsm3[H,n];


legend=Framed[SwatchLegend[{RGBColor[0.34398, 0.49112, 0.89936],RGBColor[0.97, 0.606, 0.081],RGBColor[0.91, 0.318, 0.243],RGBColor[0.448, 0.69232, 0.1538]},{"DSM0","DSM1","DSM2","DSM3"},LegendLayout->"Row"],Background->White]


hplot=Windowh[{
{(ToPureVars@dsm0dl[H,n])[[1]]},
{(ToPureVars@dsm1dl[H,n])[[1]]},
{(ToPureVars@dsm2dl[H,n])[[1]]},
{(ToPureVars@dsm3dl[H,n])[[1]]}},"",400,{0.95,1.3}]



umplot=Windowh[{{(ToPureVars@dsm0dl[H,n])[[2]]},{(ToPureVars@dsm1dl[H,n])[[2]]},{(ToPureVars@dsm2dl[H,n])[[2]]},{(ToPureVars@dsm3dl[H,n])[[2]]}},"",400,{0,2}]
a1plot=Windowh[{{ConstantArray[90,n]},{(ToPureVars@dsm1dl[H,n])[[3]]},{(ToPureVars@dsm2dl[H,n])[[3]]},{(ToPureVars@dsm3dl[H,n])[[3]]}},"",400,{-.1,.1}]
a2plot=Windowh[{{ConstantArray[90,n]},{(ToPureVars@dsm2dl[H,n])[[4]]},{(ToPureVars@dsm3dl[H,n])[[5]]}},"",400,{-.1,.1}]
wmplot=Windowh[{{(ToPureVars@dsm0dl[H,n])[[3]]},{(ToPureVars@dsm1dl[H,n])[[4]]},{(ToPureVars@dsm2dl[H,n])[[5]]},{(ToPureVars@dsm3dl[H,n])[[6]]}},"",400,{-50,50}]
g1plot=Windowh[{{ConstantArray[90,n]},{(ToPureVars@dsm1dl[H,n])[[5]]},{(ToPureVars@dsm2dl[H,n])[[6]]},{(ToPureVars@dsm3dl[H,n])[[7]]}},"",400,{-50,50}]
g2plot=Windowh[{{ConstantArray[90,n]},{(ToPureVars@dsm2dl[H,n])[[7]]},{(ToPureVars@dsm3dl[H,n])[[8]]}},"",400,{-7,7}]
qmplot=Windowh[{{(ToPureVars@dsm0dl[H,n])[[4]]},{(ToPureVars@dsm1dl[H,n])[[6]]},{(ToPureVars@dsm2dl[H,n])[[8]]},{(ToPureVars@dsm3dl[H,n])[[10]]}},"",400,{-.03,.03}]
k1plot=Windowh[{{ConstantArray[90,n]},{(ToPureVars@dsm1dl[H,n])[[7]]},{(ToPureVars@dsm2dl[H,n])[[9]]},{(ToPureVars@dsm3dl[H,n])[[11]]}},"",400,{-.015,.015}]
k2plot=Windowh[{{ConstantArray[90,n]},{(ToPureVars@dsm2dl[H,n])[[10]]},{(ToPureVars@dsm3dl[H,n])[[12]]}},"",400,{-.01,.01}]


(* export variables *)


(* hvals=Flatten[{{(ToPureVars@dsm0dl[H,n])[[1]]},{(ToPureVars@dsm1dl[H,n])[[1]]},{(ToPureVars@dsm2dl[H,n])[[1]]},{(ToPureVars@dsm3dl[H,n])[[1]]}},1];
umvals=Flatten[{{(ToPureVars@dsm0dl[H,n])[[2]]},{(ToPureVars@dsm1dl[H,n])[[2]]},{(ToPureVars@dsm2dl[H,n])[[2]]},{(ToPureVars@dsm3dl[H,n])[[2]]}},1];
a1vals=Flatten[{{(ToPureVars@dsm1dl[H,n])[[3]]},{(ToPureVars@dsm2dl[H,n])[[3]]},{(ToPureVars@dsm3dl[H,n])[[3]]}},1];
a2vals=Flatten[{{(ToPureVars@dsm2dl[H,n])[[4]]},{(ToPureVars@dsm3dl[H,n])[[5]]}},1];
wmvals=Flatten[{{(ToPureVars@dsm0dl[H,n])[[3]]},{(ToPureVars@dsm1dl[H,n])[[4]]},{(ToPureVars@dsm2dl[H,n])[[5]]},{(ToPureVars@dsm3dl[H,n])[[6]]}},1];
g1vals=Flatten[{{(ToPureVars@dsm1dl[H,n])[[5]]},{(ToPureVars@dsm2dl[H,n])[[6]]},{(ToPureVars@dsm3dl[H,n])[[7]]}},1];
g2vals=Flatten[{{(ToPureVars@dsm2dl[H,n])[[7]]},{(ToPureVars@dsm3dl[H,n])[[8]]}},1];
qmvals=Flatten[{{(ToPureVars@dsm0dl[H,n])[[4]]},{(ToPureVars@dsm1dl[H,n])[[6]]},{(ToPureVars@dsm2dl[H,n])[[8]]},{(ToPureVars@dsm3dl[H,n])[[10]]}},1];
k1vals=Flatten[{{(ToPureVars@dsm1dl[H,n])[[7]]},{(ToPureVars@dsm2dl[H,n])[[9]]},{(ToPureVars@dsm3dl[H,n])[[11]]}},1];
k2vals=Flatten[{{(ToPureVars@dsm2dl[H,n])[[10]]},{(ToPureVars@dsm3dl[H,n])[[12]]}},1]; *)


(* Export[NotebookDirectory[]<>"dsm-hierarchy/h.csv",hvals]
Export[NotebookDirectory[]<>"dsm-hierarchy/um.csv",umvals]
Export[NotebookDirectory[]<>"dsm-hierarchy/a1.csv",a1vals]
Export[NotebookDirectory[]<>"dsm-hierarchy/a2.csv",a2vals]
Export[NotebookDirectory[]<>"dsm-hierarchy/wm.csv",wmvals]
Export[NotebookDirectory[]<>"dsm-hierarchy/g1.csv",g1vals]
Export[NotebookDirectory[]<>"dsm-hierarchy/g2.csv",g2vals]
Export[NotebookDirectory[]<>"dsm-hierarchy/qm.csv",qmvals]
Export[NotebookDirectory[]<>"dsm-hierarchy/k1.csv",k1vals]
Export[NotebookDirectory[]<>"dsm-hierarchy/k2.csv",k1vals]*)



(* ::Subsubsection:: *)
(*Create Error Plot*)


L1Average=Norm[#,1]/Length[#]&;


(* h *)
hdiff=Association@MapThread[Rule,{{1,2,3},L1Average[Differences[{{(ToPureVars@dsm0dl[H,n])[[1]]},{(ToPureVars@dsm1dl[H,n])[[1]]},{(ToPureVars@dsm2dl[H,n])[[1]]},{(ToPureVars@dsm3dl[H,n])[[1]]}}][[#,1]]]/L1Average[(ToPureVars@dsm3dl[H,n])[[1]]]& /@{1,2,3}}]


(* um *)
udiff=Association@MapThread[Rule,{{1,2,3},L1Average[Differences[{{(ToPureVars@dsm0dl[H,n])[[2]]},{(ToPureVars@dsm1dl[H,n])[[2]]},{(ToPureVars@dsm2dl[H,n])[[2]]},{(ToPureVars@dsm3dl[H,n])[[2]]}}][[#,1]]]/L1Average[(ToPureVars@dsm3dl[H,n])[[2]]]& /@{1,2,3}}]


(* a1 *)
a1diff=Association@MapThread[Rule,{{2,3},L1Average[Differences[{{(ToPureVars@dsm1dl[H,n])[[3]]},{(ToPureVars@dsm2dl[H,n])[[3]]},{(ToPureVars@dsm3dl[H,n])[[3]]}}][[#,1]]]/L1Average[(ToPureVars@dsm3dl[H,n])[[3]]]& /@{1,2}}]


(* wm *)
wmdiff=Association@MapThread[Rule,{{1,2,3},L1Average[Differences[{{(ToPureVars@dsm0dl[H,n])[[3]]},{(ToPureVars@dsm1dl[H,n])[[4]]},{(ToPureVars@dsm2dl[H,n])[[5]]},{(ToPureVars@dsm3dl[H,n])[[6]]}}][[#,1]]]/L1Average[(ToPureVars@dsm3dl[H,n])[[6]]]& /@{1,2,3}}]


(* g1 *)
g1diff=Association@MapThread[Rule,{{2,3},L1Average[Differences[{{(ToPureVars@dsm1dl[H,n])[[5]]},{(ToPureVars@dsm2dl[H,n])[[6]]},{(ToPureVars@dsm3dl[H,n])[[7]]}}][[#,1]]]/L1Average[(ToPureVars@dsm3dl[H,n])[[7]]]& /@{1,2}}]


(* qm *)
qmdiff=Association@MapThread[Rule,{{1,2,3},L1Average[Differences[{{(ToPureVars@dsm0dl[H,n])[[4]]},{(ToPureVars@dsm1dl[H,n])[[6]]},{(ToPureVars@dsm2dl[H,n])[[8]]},{(ToPureVars@dsm3dl[H,n])[[10]]}}][[#,1]]]/L1Average[(ToPureVars@dsm3dl[H,n])[[10]]]& /@{1,2,3}}]


(* k1 *)
k1diff=Association@MapThread[Rule,{{2,3},L1Average[Differences[{{(ToPureVars@dsm1dl[H,n])[[7]]},{(ToPureVars@dsm2dl[H,n])[[9]]},{(ToPureVars@dsm3dl[H,n])[[11]]}}][[#,1]]]/L1Average[(ToPureVars@dsm3dl[H,n])[[11]]]& /@{1,2}}]


ListLogPlot[{hdiff,udiff,a1diff,wmdiff,g1diff,qmdiff,k1diff},Joined->True,PlotMarkers->"OpenMarkers",Frame->True,Background->White,GridLines->{{1,2,3},Automatic},
FrameLabel->{{None,"err"},{None,"N"}},AspectRatio->3/4,
PlotStyle->{Red,Purple,Blue,Cyan,Orange,Green,Brown},LabelStyle->Directive[Black,22],
(*Epilog->{Inset[Framed["nonlinear",Background->White],{2,Log[.8]}]},*)FrameTicks->{Automatic,{{1,2,3},None}},PlotRangePadding->{{.1,.25},{.25,.25}},PlotRange->{0.000001,3},
PlotLegends->{h,Subscript[u,m],Subscript[a,1],Subscript[w,m],Subscript["g",1],Subscript[q,m],Subscript[k,1]},Background->White,GridLinesStyle->Black,FrameStyle->Directive[Thick,Black]
]


(* Export the data *) 

(*Export[NotebookDirectory[]<>"diffs.csv",{hdiff,udiff,a1diff,wmdiff,g1diff,qmdiff,k1diff}]*)


(* ::Subsection:: *)
(*Linear Equations*)


(* ::Text:: *)
(*The linear equations are not giving accurate results for this experiment with the default parameters.*)


dsm0linear[H,n]=DSMSolve[varnamesLDSM0,LDSM0eqns,ConstantArray[0,Length@varnamesLDSM0]&,dictLDSM0,initSMEcompareLinear[4,h0,u0,0*-1],constrLDSM0,loadParameters,.995,tend,n];

dsm1linear[H,n]=DSMSolve[varnamesLDSM1,LDSM1eqns,ConstantArray[0,Length@varnamesLDSM1]&,dictLDSM1,initSMEcompareLinear[7,h0,u0,0*-1],
constrLDSM1,loadParameters,.995,tend,n];

dsm2linear[H,n]=DSMSolve[varnamesLDSM2,LDSM2eqns,ConstantArray[0,Length@varnamesLDSM2]&,dictLDSM2,initSMEcompareLinear[10,h0,u0,0*-1],
constrLDSM2,loadParameters,.995,tend,n];

dsm3linear[H,n]=DSMSolve[varnamesLDSM3,LDSM3eqns,ConstantArray[0,Length@varnamesLDSM3]&,dictLDSM3,initSMEcompareLinear[13,h0,u0,0*-1],
constrLDSM3,loadParameters,.995,tend,n];


dsm0lineardl[H,n]=ToDimlessVarsLinear@dsm0linear[H,n];
dsm1lineardl[H,n]=ToDimlessVarsLinear@dsm1linear[H,n];
dsm2lineardl[H,n]=ToDimlessVarsLinear@dsm2linear[H,n];
dsm3lineardl[H,n]=ToDimlessVarsLinear@dsm3linear[H,n];


(*linear plots*)


lhplot=Windowh[{
{(dsm0lineardl[H,n])[[1]]},
{(dsm1lineardl[H,n])[[1]]},
{(dsm2lineardl[H,n])[[1]]},
{(dsm3lineardl[H,n])[[1]]}},"",400,{0.95,1.3}]



lumplot=Windowh[{{(dsm0lineardl[H,n])[[2]]},{(dsm1lineardl[H,n])[[2]]},{(dsm2lineardl[H,n])[[2]]},{(dsm3lineardl[H,n])[[2]]}},"",400,{0,2}]
la1plot=Windowh[{{ConstantArray[90,n]},{(dsm1lineardl[H,n])[[3]]},{(dsm2lineardl[H,n])[[3]]},{(dsm3lineardl[H,n])[[3]]}},"",400,{-.1,.1}]
la2plot=Windowh[{{ConstantArray[90,n]},{(dsm2lineardl[H,n])[[4]]},{(dsm3lineardl[H,n])[[5]]}},"",400,{-.1,.1}]
lwmplot=Windowh[{{(dsm0lineardl[H,n])[[3]]},{(dsm1lineardl[H,n])[[4]]},{(dsm2lineardl[H,n])[[5]]},{(dsm3lineardl[H,n])[[6]]}},"",400,{-50,50}]
lg1plot=Windowh[{{ConstantArray[90,n]},{(dsm1lineardl[H,n])[[5]]},{(dsm2lineardl[H,n])[[6]]},{(dsm3lineardl[H,n])[[7]]}},"",400,{-50,50}]
lg2plot=Windowh[{{ConstantArray[90,n]},{(dsm2lineardl[H,n])[[7]]},{(dsm3lineardl[H,n])[[8]]}},"",400,{-7,7}]
qmplot=Windowh[{{(dsm0lineardl[H,n])[[4]]},{(dsm1lineardl[H,n])[[6]]},{(dsm2lineardl[H,n])[[8]]},{(dsm3lineardl[H,n])[[10]]}},"",400,{-.03,.03}]
lk1plot=Windowh[{{ConstantArray[90,n]},{(dsm1lineardl[H,n])[[7]]},{(dsm2lineardl[H,n])[[9]]},{(dsm3lineardl[H,n])[[11]]}},"",400,{-.015,.015}]
lk2plot=Windowh[{{ConstantArray[90,n]},{(dsm2lineardl[H,n])[[10]]},{(dsm3lineardl[H,n])[[12]]}},"",400,{-.01,.01}]

