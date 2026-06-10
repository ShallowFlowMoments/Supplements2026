(* ::Package:: *)

DSMSolve::usage="DSMSolve[varnames_,setEquations_,force_,vardict_,initSolution_,fdConstraints_,loadParameters_,cflno_,tEnd_,resolution_] 
solves the system with periodic BC";


(* ::Input::Initialization:: *)
SetAttributes[{P,A0},Listable];SetAttributes[poissonStep,HoldAll];
ZeroMatrix[n_]:=0IdentityMatrix[n]
Diag[list_]:=DiagonalMatrix[list]

interpolate[pts_,val_,dim_]:=Block[{len=Length[Flatten[val]]},
Interpolation[ArrayReshape[Flatten[Thread[{ArrayReshape[Flatten[pts],{len,dim}],Flatten[val]}]],{len,dim+1}],InterpolationOrder->1]
]

SetBC[list_,BCs_]:=SetRightBC[SetLeftBC[list,BCs[[1]]],BCs[[2]]];
SetRightBC[list_,bc_]:=If["i"==bc,ArrayPad[list,{0,2},"Extrapolated",InterpolationOrder->1],ArrayPad[list,{0,2},bc]];
SetLeftBC[list_,bc_]:=If["i"==bc,ArrayPad[list,{2,0},"Extrapolated",InterpolationOrder->1],ArrayPad[list,{2,0},bc]];
printStartOutput[startDate_]:=Echo["Simulation started on: "<>DateString@startDate]<>".";
printEndOutput[startDate_,endDate_,tEnd_]:=Echo["Simulation of '" <>modeltag<>"' ended on "<>DateString@endDate<> ". Simulated "<>ToString@tEnd<>" time units ("<>ToString@step<>" steps) in " <>ToString[endDate-startDate]]<>".";
L1ErrH[]:=Total[Abs[Uraw[[1]]-hinit[pts,time]]]/nx


(* ::Input::Initialization:: *)
gam=1.5;
Rlimit[dm_,dp_]=If[dp dm<=0,0,Sign[dp]Min[Abs[gam dm],Abs[(2dp+dm)/3],Abs[(1+gam)/2 dp]]]; (*Limiter*)
Rfull[dm_,dp_]=Sign[dp]Max[0,Min[Sign[dp] (2dp+dm)/3,Max[-Sign[dp]dm,Min[gam Sign[dp]dm,Sign[dp] (2dp+dm)/3,(1+gam)/2 Abs[dp]]]]];

(*Recon[list_,opt_]:=Map[Map[Function[{um,u,up},{u+1/2Rlimit[u-um,up-u],u-1/2Rlimit[up-u,u-um]}]@@#&,Partition[SetBC[#,{"i","i"}],3,1]]&,list]*)
Recon[list_,"Fixed"]:=MapThread[Map[Function[{um,u,up},{u+1/2 Rlimit[u-um,up-u],u-1/2 Rlimit[up-u,u-um]}]@@#&,Partition[SetBC[#1,#2],3,1]]&,{list,BClist}]
Recon[list_,"Periodic"]:=Map[Map[Function[{um,u,up},{u+1/2 Rlimit[u-um,up-u],u-1/2 Rlimit[up-u,u-um]}]@@#&,Partition[ArrayPad[#,2,"Periodic",InterpolationOrder->1],3,1]]&,list]


(* ::Input::Initialization:: *)
Residuum[U_]:=Block[{Urecons1,Uright1,Uleft1,CMax,delta=0.95},
Urecons=Recon[U,bcType];

Uright=Urecons[[All,;;-2,1]];
Uleft=Urecons[[All,2;;,2]];
dU=1/2 Differences[ArrayPad[U,{{0,0},{1,1}},bcType],{0,1},2];
CMax=ConstantArray[ListConvolve[{1/2,1/2},MaxSpeed@@U,{1,-1},"Periodic"],Length[U]];
cMax=Max[CMax];1/dx(-Differences[1/2 (Flux@@Uright+Flux@@Uleft)+delta/2 cMax(Uright-Uleft),{0,1}]+ Transpose[MapThread[#1 . #2&,{A0@@U,Transpose[dU]}]]+Transpose[dT@@#&/@Transpose[U] dxTimesdHbList])+Transpose[P@@U]+Transpose[T@@#&/@Transpose[U] HbList]+(Transpose@MapThread[(#1 . #2)&,{F@@@(Transpose@Uraw),forceList[tmin]}])
]


(* ::Input::Initialization:: *)
FiniteVolumeRun[nx_,tend_]:=Block[{H1,HU1,H2,HU2},
startDate=DateObject[];
printStartOutput[startDate];


tmin=0;


forceList[t_]=force[#,t]& /@pts;
(*forceList[t_]=force[#,t]& /@pts;*)
Uraw=Transpose[Table[init[xs[i],time],{i,1,nx}]];
U1=Uraw;
U2=Uraw;
(*qvals=-Uraw[[1]]-HbList; (*initialize with -(h+hb) *)
(*hInterfaceVals=ListConvolve[{-1/16,9/16,9/16,-1/16},Uraw[[1]],-2];*)
pvals =LinearSolve[Normal@A/.{hb[i_]\[Rule]Indexed[HbList,i],ddhb[i_]\[Rule]Indexed[ddHbList,i],h[i_]\[Rule]Indexed[Uraw[[1]],i],t\[Rule]time},-Normal@b/.{hb[i_]\[Rule]Indexed[HbList,i],ddhb[i_]\[Rule]Indexed[ddHbList,i],h[i_]\[Rule]Indexed[Uraw[[1]],i],t\[Rule]time}];
(*pvals=ListConvolve[{-1/16,9/16,9/16,-1/16},pInterfaceVals,2];*)
dxpvals=ListConvolve[{-1/12,8/12,0,-8/12,1/12},pvals,3]/(dx);
qvals=qvals+pvals;
Uraw[[4]]=qvals;*)
poissonStep[Uraw];

cMax=0.0;
dt=0.1dx;
If[time+dt>=tend&&time<tend,dt=tend-time+10^-8];

step=0;
Uenlarged=Join[Uraw,{HbList,dHbList}];
times={0};
(*perturbations={Total/@Abs[Uraw-stationary[pts]]};times={0};*)
(*ePot=Epot@@Uenlarged;
eKin=Ekin@@Uenlarged;
Ein=u0(eKin[[1]]+ePot[[1]]);
Eout=u0(eKin[[nx]]+ePot[[nx]]);
puin=h0 pu@@Uenlarged[[;;,1]];
puout=h0 pu@@Uenlarged[[;;,nx]];
eSource=res@@Uenlarged;

dtE=Eout-Ein+puout-puin-Total[eSource]/nx;
eSourceList={Total[eSource]/nx};
energies={{Total[ePot+eKin]/nx,0,Total[ePot+eKin]/nx,Ein,puin,Eout,puout,dtE}};*)
While[time<tend,

(*stage 1*)
tmin=time;



U1[[1;;(nVars-nPVars)]]=(Uraw+dt Residuum[Uraw])[[1;;(nVars-nPVars)]];

(*stage 2*)
tmin=time+dt;

(*u1pvals =LinearSolve[Normal@A/.{hb[i_]\[Rule]Indexed[HbList,i],ddhb[i_]\[Rule]Indexed[ddHbList,i],h[i_]\[Rule]Indexed[U1[[1]],i],t\[Rule]tmin},-Normal@b/.{hb[i_]\[Rule]Indexed[HbList,i],ddhb[i_]\[Rule]Indexed[ddHbList,i],h[i_]\[Rule]Indexed[U1[[1]],i],t\[Rule]tmin}];
(*u1pvals=ListConvolve[{-1/16,9/16,9/16,-1/16},u1pInterfaceVals,2];*)
u1dxpvals=ListConvolve[{-1/12,8/12,0,-8/12,1/12},u1pvals,3]/(dx);
pvals=u1pvals;
dxpvals=u1dxpvals;

qvals=pvals-U1[[1]]-HbList;
U1[[4]]=qvals;*)
poissonStep[U1];

U2[[1;;(nVars-nPVars)]]=(3/4 Uraw+1/4 (U1+dt Residuum[U1]))[[1;;(nVars-nPVars)]];

(*stage 3*)
tmin=time+1/2 dt;

(*u2pvals =LinearSolve[Normal@A/.{hb[i_]\[Rule]Indexed[HbList,i],ddhb[i_]\[Rule]Indexed[ddHbList,i],h[i_]\[Rule]Indexed[U2[[1]],i],t\[Rule]tmin},-Normal@b/.{hb[i_]\[Rule]Indexed[HbList,i],ddhb[i_]\[Rule]Indexed[ddHbList,i],h[i_]\[Rule]Indexed[U2[[1]],i],t\[Rule]tmin}];
(*u2pvals=ListConvolve[{-1/16,9/16,9/16,-1/16},u2pInterfaceVals,2];*)
u2dxpvals=ListConvolve[{-1/12,8/12,0,-8/12,1/12},u2pvals,3]/(dx);
pvals=u2pvals;
dxpvals=u2dxpvals;

qvals=pvals-U2[[1]]-HbList;
U2[[4]]=qvals ;*)
poissonStep[U2];


Uraw[[1;;(nVars-nPVars)]]=(1/3 Uraw+2/3 (U2+dt Residuum[U2]))[[1;;(nVars-nPVars)]];

tmin=time+dt;
(*
urawpvals =LinearSolve[Normal@A/.{hb[i_]\[Rule]Indexed[HbList,i],ddhb[i_]\[Rule]Indexed[ddHbList,i],h[i_]\[Rule]Indexed[Uraw[[1]],i],t\[Rule]tmin},-Normal@b/.{hb[i_]\[Rule]Indexed[HbList,i],ddhb[i_]\[Rule]Indexed[ddHbList,i],h[i_]\[Rule]Indexed[Uraw[[1]],i],t\[Rule]tmin}];
(*urawpvals=ListConvolve[{-1/16,9/16,9/16,-1/16},urawpInterfaceVals,2];*)
urawdxpvals=ListConvolve[{-1/12,8/12,0,-8/12,1/12},urawpvals,3]/(dx);
pvals=urawpvals;
dxpvals=urawdxpvals;

qvals=pvals-Uraw[[1]]-HbList;
Uraw[[4]]=qvals ;*)
poissonStep[Uraw];

Uenlarged=Join[Uraw,{HbList,dHbList}];

(*ePot=Epot@@Uenlarged;
eKin=Ekin@@Uenlarged;
Ein=u0(eKin[[1]]+ePot[[1]]);
Eout=u0(eKin[[nx]]+ePot[[nx]]);
puin=h0 pu@@Uenlarged[[;;,1]];
puout=h0 pu@@Uenlarged[[;;,nx]];
eSource=res@@Uenlarged;
dtE=Eout-Ein+puout-puin-Total[eSource]/nx;

If[Mod[step,resUpdateInterval]\[Equal]0,{
eSourceList=Append[eSourceList,Total[eSource]/nx];
energies=Append[energies,{Total[ePot+eKin]/nx,Total[eSourceList]*resUpdateInterval dt,Total[ePot+eKin]/nx-Total[eSourceList]*resUpdateInterval dt,Ein,puin,Eout,puout,dtE,Total[eSource]/nx}];
(*perturbations=Append[perturbations,Total/@Abs[Uraw-stationary[pts]]/nx];*)
times=Append[times,time];}];*)
If[Min[h0+Uraw[[1]]]<=0,Echo["WARNING: Total height below zero!!"]];

time+=dt;
step++;
Capture[];

dt=CFL dx/cMax;
(*CFL=cMax dt/dx;*)
If[time+dt>=tend&&time<tend,dt=tend-time+10^-8];
];
endDate=DateObject[];
printEndOutput[startDate,endDate,tend];
{time,step}
];




(* ::Input::Initialization:: *)
initFDsolution=Function[{fdConstraints},

varsFull=Flatten[Table[{pvars[[j]][i]},{i,1,nx},{j,1,nPVars}]];
(*xs[i_]=x1+(i-0.5)dx;*)
eqns=Flatten[fdConstraints[nx]]/. {a_[nx+1]/;MemberQ[{p,h,q,w,u,k1,k2,k3,a1,a2,a3,g1,g2,g3,hb,dhb,ddhb},a]->a[1],
a_[0]/;MemberQ[{p,h,q,w,u,k1,k2,k3,a1,a2,a3,g1,g2,g3,hb,dhb,ddhb},a]->a[nx],
a_[nx+2]/;MemberQ[{p,h,q,w,u,k1,k2,k3,a1,a2,a3,g1,g2,g3,hb,dhb,ddhb},a]->a[2],
a_[-1]/;MemberQ[{p,h,q,w,u,k1,k2,k3,a1,a2,a3,g1,g2,g3,hb,dhb,ddhb},a]->a[nx-1]}(* /. {hb[i_]->Indexed[HbList,i],
dhb[i_]->Indexed[dHbList,i],
ddhb[i_]->Indexed[ddHbList,i]}*);
{b,A}=CoefficientArrays[eqns,varsFull];
bhat=-b;
];


(* ::Input::Initialization:: *)
basedict[Uvec_]:={hb[i_]->Indexed[HbList,i],
dhb[i_]->Indexed[dHbList,i],
ddhb[i_]->Indexed[ddHbList,i],
h[i_]->Indexed[Uvec[[1]],i],
t->tmin};

poissonStep[U_]:=
If[constrIsEmpty,Return[],(
	ucurrentsol =LinearSolve[Normal@A/.basedict[U]/.dict[U],bhat/.basedict[U]/.dict[U](*,Method->{"Krylov","Method"->"BiCGSTAB","Preconditioner"->"ILU0","Tolerance"->10^-6}*)] ;
upressurevals=(ucurrentsol[[#;;nPVars nx-(nPVars-#);;nPVars]])&/@Range[nPVars];
(*u1pvals=ListConvolve[{-1/16,9/16,9/16,-1/16},u1pInterfaceVals,2];*)
(*u1dxpvals=ListConvolve[{-1/12,8/12,0,-8/12,1/12},u1pvals,3]/(dx);*)
(*pvals=upvals;*)
(*dxpvals=u1dxpvals;*)
kappavals=upressurevals[[1;;nPVars-1]];
uinsert=Join[{upressurevals[[nPVars]]-(U[[1]]+HbList+Total[kappavals])},kappavals];
(*transform p into q*)
(*Echo[ucurrentsol];*)

(*idea: shift depending on last entry of varnames *)
U[[(nVars)-nPVars+1;;(nVars)]]=uinsert;
(*qvals=pvals-U[[1]]-HbList;
U[[4]]=qvals;*)
);];


(* ::Input::Initialization:: *)
DSMSolve[varnames_,setEquations_,forcearg_,vardict_,initSolution_,fdConstraints_,loadParameters_,cflno_,tEnd_,resolution_] :=(
Clear[Uraw];
g=1;nx=resolution;
time=0;
dx=(x2-x1)/nx;
xs[i_]=x1+(i-0.5)dx;
pts=Table[xs[i],{i,1,nx}];


HbList=Hb /@pts;
dxTimesdHbList=1/2 Differences[ArrayPad[{HbList},{{0,0},{1,1}},"Periodic"],{0,1},2][[1]];
dHbList=1/dx Table[-(-1/(12)hb[i-2]+8/(12)hb[i-1]-8/(12)hb[i+1]+1/(12)hb[i+2]),{i,1,nx}]/. {hb[nx+1]->hb[1],hb[0]->hb[nx],hb[nx+2]->hb[2],hb[-1]->hb[nx-1]} /. hb[i_]->Indexed[HbList,i];
ddHbList=Table[(-hb[i+2]+16hb[i+1]-30hb[i]+16hb[i-1]-hb[i-2])/(12 dx^2),{i,1,nx}]/. {hb[nx+1]->hb[1],hb[0]->hb[nx],hb[nx+2]->hb[2],hb[-1]->hb[nx-1]} /. hb[i_]->Indexed[HbList,i];

force=forcearg;
ref=initSolution;
nVars=If[varnames[[-1]]=="h_b",Length[varnames]-1,Length[varnames]];

constrIsEmpty=ListQ[fdConstraints]&&{}==fdConstraints;
nPVars=If[2==nVars,0,Quotient[nVars,3]];
pvars=If[2==nVars,{},Join[Symbol/@varnames[[(nVars)-nPVars+2;;(nVars)]],{p}]];
dict=vardict;
(*loadParameters[];(*loadInitialData[];*)*)
       (*Fr=u0/Sqrt[h0];*)Echo["Fr="<>ToString@Fr];
CFL=cflno;
(*S=h0; !!!!!! what about factor 2 because lengh of domain is 2??? *) 
       modeltag="DSMSolve";
 (*the larger the value, the slower the simulation*)
setEquations[];
(*Flux[h_,u_,w_,q_]={u0 h+h0 u,u0 u+h+q,u0 w,0};
P[h_,u_,w_,q_]={0,0,1/h0 q, 0}; (*source term from flow vars*)
T[h_,u_,w_,q_]={0 ,0,0,0}; (*source term from hb*)
dT[h_,u_,w_,q_]={0 ,-g,0,0}; (*source term from dx hb*)
A0[h_,u_,w_,q_]=Block[{},0IdentityMatrix[4]];*)

init=initSolution;

If[!constrIsEmpty,initFDsolution[fdConstraints]];

(*Epot[h_,u_,w_,q_,hb_]=1/2+h+1/2h^2;
Ekin[h_,u_,w_,q_]= (1/2 (u+1)^2+1/2w^2)+h(1/2 (u+1)^2+1/2w^2);*)
(*Epot[h_,u_,w_,q_,hb_,dxhb_]=1/2 h^2 + h hb;
       Ekin[h_,u_,w_,q_,hb_,dxhb_]=h0 (1/2 u^2 + 1/2 w^2);
pu[h_,u_,w_,q_,hb_,dxhb_]= (h+hb)u+q u;
res[h_,u_,w_,q_,hb_,dxhb_]=u0(h+q)dxhb;*)
(*:={{0,"i"},{"i",0},{0,"i"},{0,"i"}};*)
BClist:={{0,"i"},{"i",eps((Uraw[[1,nx]])+qvals[[nx]])},{0,"i"}};
(*BClist:={{0,"i"},{"i",eps Uraw[[4,nx]]},{0,"i"},{ -eps (Uraw[[2,1]]),"i"}};*)
(*BClist:={'
{Uraw[[4,1]]/cr^2,"i"},
{"i",(g Uraw[[1,nx]] +Uraw[[4,nx]])/Sqrt[(cr^2+g)h0]},
{0,"i"},
{-g Uraw[[1,1]] - Sqrt[(cr^2+g)h0] Uraw[[2,1]],"i"}};*)
(*varnames={"h","u_m","w_m","q_m"};*)
(*stationary[x_]={(h[x]), (um[x]),wm[x]}/.dsm0;
init[x_]=Join[{stationary[x][[1]]+ h0 mu (1+Cos[\[Pi] x])(*Sin[\[Pi] x]*)},stationary[x][[2;;3]]];*)





FiniteVolumeRun[nx,tEnd];
Uraw
);

(*MakeGif[];FinishCapture[];*)


