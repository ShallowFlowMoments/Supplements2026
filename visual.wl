(* ::Package:: *)

(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*Visualization for Dynamic Output*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Input:: *)
(*iSize=400;*)
(*color={Red,Green,Blue,Purple};*)
(*markers={\[FilledUpTriangle],\[FilledSquare],\[FilledDiamond],\[FilledCircle],\[FivePointedStar]};*)
(*padding={{40,10},{20,20}}*)


(* ::Input::Initialization:: *)
VisualList[UList_,dispnames_,iSize_,ranges_:{{0,1},{0,1},{-.5,.5},{-.5,.5}}]:=Block[{ny},
Grid[{{
Show[Plot[{Hb[x]+ref[x,time][[1]],Hb[x]},{x,x1,x2},LabelStyle->Directive[Black, FontSize->22],FrameStyle->
Directive[Black,Thick],ImageSize->Medium,PlotLabel->dispnames[[1]]<>"(x,t)"(* at time: "<>ToString[time]<>" (step: "<>ToString[step]<>", \[CapitalDelta]t = "<>ToString[dt]<>", CFL = "<>ToString[CFL]<>")"*),PlotRange->{{x1,x2},ranges[[1]]},ImagePadding->{{100,20},{30,20}},Frame->True,ImageSize->iSize],
MapThread[ListPlot[{pts,#1[[1]]}\[Transpose],PlotStyle->ColorData[106,#2],ImageSize->iSize]& ,{UList,Range[Length@UList]}]],
Show[Plot[ref[x,time][[2]]/ref[x,time][[1]],{x,x1,x2},LabelStyle->Directive[Black, FontSize->22],FrameStyle->
Directive[Black,Thick],ImageSize->Medium,PlotStyle->{Red,Blue,Green},PlotLabel->dispnames[[2]]<>"_Ave(x,t)",PlotRange->{{x1,x2},ranges[[2]]},ImagePadding->{{100,20},{30,20}},Frame->True,ImageSize->iSize],
MapThread[ListPlot[{pts,#1[[2]]}\[Transpose],PlotStyle->ColorData[106,#2],ImageSize->iSize]& ,{UList,Range[Length@UList]}]]},
If[Length@dispnames>2,
{Show[Plot[ref[x,time][[3]]/ref[x,time][[1]],{x,x1,x2},LabelStyle->Directive[Black, FontSize->22],FrameStyle->
Directive[Black,Thick],ImageSize->Medium,PlotLabel->dispnames[[3]]<>"_Ave(x,t) ",PlotRange->{{x1,x2},ranges[[3]]},ImagePadding->{{100,20},{30,20}},Frame->True,ImageSize->iSize],
MapThread[ListPlot[{pts,#1[[3]]}\[Transpose],PlotStyle->ColorData[106,#2],ImageSize->iSize]& ,{UList,Range[Length@UList]}]],
Show[Plot[ref[x,time][[4]],{x,x1,x2},LabelStyle->Directive[Black, FontSize->22],FrameStyle->
Directive[Black,Thick],ImageSize->Medium,PlotStyle->{Red,Blue,Green},PlotLabel->dispnames[[4]]<>"_Ave(x,t)",PlotRange->{{x1,x2},ranges[[4]]},
ImagePadding->{{100,20},{30,20}},Frame->True,ImageSize->iSize],
MapThread[ListPlot[{pts,#1[[4]]}\[Transpose],PlotStyle->ColorData[106,#2],ImageSize->iSize]& ,{UList,Range[Length@UList]}]]},{None,None}]
(*,{Show[Plot[0,{x,x1,x2},PlotLabel\[Rule]"Epot(x,t)",PlotRange\[Rule]{{x1,x2},{-.2,.5}},Frame\[Rule]True,ImageSize\[Rule]iSize,Epilog\[Rule]{Inset[Panel[modeltag],Scaled[{.9,.9}]]}],
ListPlot[Transpose[{pts,ePot}],ImageSize\[Rule]iSize]],
Show[Plot[{0},{x,x1,x2},PlotStyle\[Rule]{Red,Blue,Green},PlotLabel\[Rule]"Ekin(x,t)",PlotRange\[Rule]{{x1,x2},{-.2,1}},Frame\[Rule]True,ImageSize\[Rule]iSize,Epilog\[Rule]{Inset[Panel[modeltag],Scaled[{.9,.9}]]}],
ListPlot[Transpose[{pts,eKin}],ImageSize\[Rule]iSize]]}*)}]

];

Visual2[U_,iSize_]:=Block[{ny},
Grid[{(*{ListPlot[TemporalData[perturbations[[1;;Length@times,All]],{times}],ImageSize\[Rule]{2 iSize,267},AspectRatio\[Rule]267/800,PlotLabel\[Rule]"L1-diff from stationary; nx="<>ToString@nx<>";eps="<>ToString@eps<>";mu="<>ToString@mu<>";\[Eta]="<>ToString@\[Eta]<>";\[Sigma]="<>ToString@\[Sigma]<>";Fr="<>ToString@Fr<>";S="<>ToString@S<>";cr="<>ToString@cr,PlotLegends\[Rule]Placed[LineLegend[varnames[[All]]],Right],AxesLabel\[Rule]{"t (seconds)","L1"},LabelStyle\[Rule]"Text",Epilog\[Rule]{Inset[Panel[modeltag],Scaled[{.9,.9}]]},GridLines\[Rule]Automatic
]None,ListPlot[TemporalData[energies[[1;;Length@times,1;;3]],{times}],ImageSize\[Rule]{2 iSize,267},AspectRatio\[Rule]267/800,PlotLabel\[Rule]"Total energy; nx="<>ToString@nx<>";eps="<>ToString@eps<>";mu="<>ToString@mu<>";\[Eta]="<>ToString@\[Eta]<>";\[Sigma]="<>ToString@\[Sigma]<>";Fr="<>ToString@Fr<>";S="<>ToString@S<>";cr="<>ToString@cr,PlotLegends\[Rule]Placed[LineLegend[{"eKin+ePot","eSrc","e-eSrc"}],Right],AxesLabel\[Rule]{"t (seconds)","E"},LabelStyle\[Rule]"Text",Epilog\[Rule]{Inset[Panel[modeltag],Scaled[{.9,.9}]]},GridLines\[Rule]Automatic
]},{
ListPlot[TemporalData[energies[[1;;Length@times,4;;5]],{times}],ImageSize\[Rule]{2 iSize,267},AspectRatio\[Rule]267/800,PlotLabel\[Rule]"Inflow; nx="<>ToString@nx<>";eps="<>ToString@eps<>";mu="<>ToString@mu<>";\[Eta]="<>ToString@\[Eta]<>";\[Sigma]="<>ToString@\[Sigma]<>";Fr="<>ToString@Fr<>";S="<>ToString@S<>";cr="<>ToString@cr,PlotLegends\[Rule]Placed[LineLegend[{"Ein","puin"}],Right],AxesLabel\[Rule]{"t (seconds)","E"},LabelStyle\[Rule]"Text",Epilog\[Rule]{Inset[Panel[modeltag],Scaled[{.9,.9}]]},GridLines\[Rule]Automatic],
ListPlot[TemporalData[energies[[1;;Length@times,6;;7]],{times}],ImageSize\[Rule]{2 iSize,267},AspectRatio\[Rule]267/800,PlotLabel\[Rule]"Outflow; nx="<>ToString@nx<>";eps="<>ToString@eps<>";mu="<>ToString@mu<>";\[Eta]="<>ToString@\[Eta]<>";\[Sigma]="<>ToString@\[Sigma]<>";Fr="<>ToString@Fr<>";S="<>ToString@S<>";cr="<>ToString@cr,PlotLegends\[Rule]Placed[LineLegend[{"Eout","puOut"}],Right],AxesLabel\[Rule]{"t (seconds)","E"},LabelStyle\[Rule]"Text",Epilog\[Rule]{Inset[Panel[modeltag],Scaled[{.9,.9}]]},GridLines\[Rule]Automatic]
},{ListPlot[TemporalData[energies[[1;;Length@times,8]],{times}],ImageSize\[Rule]{2 iSize,267},AspectRatio\[Rule]267/800,PlotLabel\[Rule]"-dtE ; nx="<>ToString@nx<>";eps="<>ToString@eps<>";mu="<>ToString@mu<>";\[Eta]="<>ToString@\[Eta]<>";\[Sigma]="<>ToString@\[Sigma]<>";Fr="<>ToString@Fr<>";S="<>ToString@S<>";cr="<>ToString@cr,PlotLegends\[Rule]Placed[LineLegend[{"Flow over boundary"}],Right],AxesLabel\[Rule]{"t (seconds)","E"},LabelStyle\[Rule]"Text",Epilog\[Rule]{Inset[Panel[modeltag],Scaled[{.9,.9}]]},GridLines\[Rule]Automatic]}*)}]

];


(* ::Input::Initialization:: *)
(*Export Visualization to GIF*)
Capture[]:=If[capture,AppendTo[frames,Visual[Uraw,400]]]
MakeGif[]:=Export[NotebookDirectory[]<>modeltag<>DateString["ISODateTime"]<>".gif",ArrayPad[frames,20,"Fixed"]];
InitCapture[]:=(capture=True;frames={});
FinishCapture[]:=(capture=False;Clear@frames);


Visual[U_,iSize_]:=Block[{ny},

Grid[{{
Show[Plot[{ref[x,time][[1]],Hb[x]},{x,x1,x2},PlotLabel->"h(x,t) at time: "<>ToString[time]<>" (step: "<>ToString[step]<>", \[CapitalDelta]t = "<>ToString[dt]<>", CFL = "<>ToString[CFL]<>")",PlotRange->{{x1,x2},{.8,1.2}},Frame->True,ImageSize->iSize,Epilog->{Inset[Panel[modeltag],Scaled[{.9,.9}]]}],
ListPlot[Transpose[{pts,HbList/h0+U[[1]]}],ImageSize->iSize]],
Show[Plot[ref[x,time][[2]],{x,x1,x2},PlotStyle->{Red,Blue,Green},PlotLabel->"u_Ave(x,t) or hu_Ave(x,t)",PlotRange->{{x1,x2},{0.6,1.6}},Frame->True,ImageSize->iSize,Epilog->{Inset[Panel[modeltag],Scaled[{.9,.9}]]}],
ListPlot[Transpose[{pts,U[[2]]}],ImageSize->iSize]]},
If[Length[varnames]>2,{Show[Plot[winit[x,time],{x,x1,x2},PlotLabel->"w_Ave(x,t) ",PlotRange->{{x1,x2},{-2,2}},Frame->True,ImageSize->iSize,Epilog->{Inset[Panel[modeltag],Scaled[{.9,.9}]]}],
ListPlot[Transpose[{pts,U[[4]]}],ImageSize->iSize]],
Show[Plot[qinit[x,time],{x,x1,x2},PlotStyle->{Red,Blue,Green},PlotLabel->"q_Ave(x,t)",PlotRange->{{x1,x2},{-.4,.4}},Frame->True,ImageSize->iSize,Epilog->{Inset[Panel[modeltag],Scaled[{.9,.9}]]}],
ListPlot[Transpose[{pts,U[[6]]}],ImageSize->iSize]]},{None,None}]
(*,{Show[Plot[0,{x,x1,x2},PlotLabel\[Rule]"Epot(x,t)",PlotRange\[Rule]{{x1,x2},{-.2,.5}},Frame\[Rule]True,ImageSize\[Rule]iSize,Epilog\[Rule]{Inset[Panel[modeltag],Scaled[{.9,.9}]]}],
ListPlot[Transpose[{pts,ePot}],ImageSize\[Rule]iSize]],
Show[Plot[{0},{x,x1,x2},PlotStyle\[Rule]{Red,Blue,Green},PlotLabel\[Rule]"Ekin(x,t)",PlotRange\[Rule]{{x1,x2},{-.2,1}},Frame\[Rule]True,ImageSize\[Rule]iSize,Epilog\[Rule]{Inset[Panel[modeltag],Scaled[{.9,.9}]]}],
ListPlot[Transpose[{pts,eKin}],ImageSize\[Rule]iSize]]}*)}]

];


Window[UList_,dispname_,iSize_,range_]:=Show[Plot[{},{x,x1,x2},LabelStyle->Directive[Black, FontSize->22],FrameStyle->
Directive[Black,Thick],ImageSize->Medium,PlotLabel->dispname<>"(x,t)"(* at time: "<>ToString[time]<>" (step: "<>ToString[step]<>", \[CapitalDelta]t = "<>ToString[dt]<>", CFL = "<>ToString[CFL]<>")"*),PlotRange->{{x1,x2},range},ImagePadding->{{100,20},{30,20}},Frame->True,ImageSize->iSize],
MapThread[ListPlot[{pts,#1[[1]]}\[Transpose],PlotStyle->ColorData[106,#2],ImageSize->iSize]& ,{UList,Range[Length@UList]}]]


WindowWithRef[UList_,dispname_,iSize_,range_,plotfn_]:=Show[Plot[plotfn,{x,x1,x2},LabelStyle->Directive[Black, FontSize->22],FrameStyle->
Directive[Black,Thick],ImageSize->Medium,PlotLabel->dispname<>"(x,t)"(* at time: "<>ToString[time]<>" (step: "<>ToString[step]<>", \[CapitalDelta]t = "<>ToString[dt]<>", CFL = "<>ToString[CFL]<>")"*),PlotRange->{{x1,x2},range},ImagePadding->{{100,20},{30,20}},Frame->True,ImageSize->iSize],
MapThread[ListPlot[{pts,#1[[1]]}\[Transpose],PlotStyle->ColorData[106,#2],ImageSize->iSize]& ,{UList,Range[Length@UList]}]]
