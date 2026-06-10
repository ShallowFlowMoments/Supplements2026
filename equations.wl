(* ::Package:: *)

LDSM0eqns=Function[{},(
MaxSpeed[h_,u_,w_,q_]:=ConstantArray[Max[Abs[-Sqrt[g h0]+u0],Abs[Sqrt[g h0]+u0]],nx];
Flux[h_,u_,w_,q_]={u0 h+h0 u,u0 u+h+q,u0 w,0};
P[h_,u_,w_,q_]={0,0,1/h0 q, 0}; (*source term from flow vars*)
T[h_,u_,w_,q_]={0 ,0,0,0}; (*source term from hb*)
dT[h_,u_,w_,q_]={0 ,-g,0,0}; (*source term from dx hb*)
A0[h_,u_,w_,q_]=Block[{},0IdentityMatrix[4]];
F[h_,u_,w_,q_]=DiagonalMatrix[{1,1,1,0}];
);];


LDSM1eqns=Function[{},(
MaxSpeed[h_,u_,\[Alpha]1_,w_,\[Gamma]1_,q_,\[Kappa]1_]:=ConstantArray[Max[Abs[-Sqrt[g h0]+u0],Abs[Sqrt[g h0]+u0]],nx];
Flux[h_,u_,\[Alpha]1_,w_,\[Gamma]1_,q_,\[Kappa]1_]={u0 h+h0 u,u0 u+h+q,u0 \[Alpha]1+\[Kappa]1,u0 w,u0 \[Gamma]1,0,0};
P[h_,u_,\[Alpha]1_,w_,\[Gamma]1_,q_,\[Kappa]1_]={0,0,0,1/h0(q+\[Kappa]1),3/h0 (-q+\[Kappa]1),0,0}; (*source term from flow vars*)
T[h_,u_,\[Alpha]1_,w_,\[Gamma]1_,q_,\[Kappa]1_]={0 ,0,0,0,0,0,0}; (*source term from hb*)
dT[h_,u_,\[Alpha]1_,w_,\[Gamma]1_,q_,\[Kappa]1_]={0,-g,0,0,0,0,0}; (*source term from dx hb*)
A0[h_,u_,\[Alpha]1_,w_,\[Gamma]1_,q_,\[Kappa]1_]=Block[{},0IdentityMatrix[7]];
F[h_,u_,\[Alpha]1_,w_,\[Gamma]1_,q_,\[Kappa]1_]=DiagonalMatrix[{1,1,1,1,1,0,0}];
);];


LDSM2eqns=Function[{},(
MaxSpeed[h_,u_,\[Alpha]1_,\[Alpha]2_,w_,\[Gamma]1_,\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]:=ConstantArray[Max[Abs[-Sqrt[g h0]+u0],Abs[Sqrt[g h0]+u0]],nx];
Flux[h_,u_,\[Alpha]1_,\[Alpha]2_,w_,\[Gamma]1_,\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]={u0 h+h0 u,u0 u+h+q,u0 \[Alpha]1+\[Kappa]1,u0 \[Alpha]2+\[Kappa]2,u0 w,u0 \[Gamma]1,u0 \[Gamma]2,0,0,0};
P[h_,u_,\[Alpha]1_,\[Alpha]2_,w_,\[Gamma]1_,\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]={0,0,0,0,1/h0(q+\[Kappa]1+\[Kappa]2),3/h0 (-q+\[Kappa]1+\[Kappa]2),5/h0(q-\[Kappa]1+\[Kappa]2),0,0,0}; (*source term from flow vars*)
T[h_,u_,\[Alpha]1_,\[Alpha]2_,w_,\[Gamma]1_,\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]=ConstantArray[0,10]; (*source term from hb*)
dT[h_,u_,\[Alpha]1_,\[Alpha]2_,w_,\[Gamma]1_,\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]={0,-g,0,0,0,0,0,0,0,0}; (*source term from dx hb*)
A0[h_,u_,\[Alpha]1_,\[Alpha]2_,w_,\[Gamma]1_,\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]=Block[{},0IdentityMatrix[10]];
F[h_,u_,\[Alpha]1_,\[Alpha]2_,w_,\[Gamma]1_,\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]=DiagonalMatrix[{1,1,1,1,1,1,1,0,0,0}];
);];


LDSM3eqns=Function[{},(
MaxSpeed[h_,u_,\[Alpha]1_,\[Alpha]2_,\[Alpha]3_,w_,\[Gamma]1_,\[Gamma]2_,\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]:=ConstantArray[Max[Abs[-Sqrt[g h0]+u0],Abs[Sqrt[g h0]+u0]],nx];
Flux[h_,u_,\[Alpha]1_,\[Alpha]2_,\[Alpha]3_,w_,\[Gamma]1_,\[Gamma]2_,\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]={u0 h+h0 u,u0 u+h+q,u0 \[Alpha]1+\[Kappa]1,u0 \[Alpha]2+\[Kappa]2,u0 \[Alpha]3+\[Kappa]3,u0 w,u0 \[Gamma]1,u0 \[Gamma]2,u0 \[Gamma]3,0,0,0,0};
P[h_,u_,\[Alpha]1_,\[Alpha]2_,\[Alpha]3_,w_,\[Gamma]1_,\[Gamma]2_,\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]={0,0,0,0,0,1/h0(q+\[Kappa]1+\[Kappa]2+\[Kappa]3),3/h0 (-q+\[Kappa]1+\[Kappa]2+\[Kappa]3),5/h0(q-\[Kappa]1+\[Kappa]2+\[Kappa]3),7/h0(-q+\[Kappa]1-\[Kappa]2+\[Kappa]3),0,0,0,0}; (*source term from flow vars*)
T[h_,u_,\[Alpha]1_,\[Alpha]2_,\[Alpha]3_,w_,\[Gamma]1_,\[Gamma]2_,\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]=ConstantArray[0,13]; (*source term from hb*)
dT[h_,u_,\[Alpha]1_,\[Alpha]2_,\[Alpha]3_,w_,\[Gamma]1_,\[Gamma]2_,\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]={0,-g,0,0,0,0,0,0,0,0,0,0,0}; (*source term from dx hb*)
A0[h_,u_,\[Alpha]1_,\[Alpha]2_,\[Alpha]3_,w_,\[Gamma]1_,\[Gamma]2_,\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]=Block[{},0IdentityMatrix[13]];
F[h_,u_,\[Alpha]1_,\[Alpha]2_,\[Alpha]3_,w_,\[Gamma]1_,\[Gamma]2_,\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]=DiagonalMatrix[{1,1,1,1,1,1,1,1,1,0,0,0,0}];
);];


DLNLSWeqns=Function[{},(
MaxSpeed[h_,hu_]:=Abs[hu/h]+1/Fr Sqrt[g h];
Flux[h_,hu_]={hu,hu^2/h+1/Fr^2 h^2/2};
P[h_,hu_]={0,0}; (*source term from flow vars*)
T[h_,hu_]={0,0}; (*source term from hb*)
dT[h_,hu_]={0,-h g}; (*!!!!!! not scaled yet !!!!!!! source term from dx hb*)
A0[h_,hu_]={{0,0},{0,0}};
F[h_,hu_]={{1,0},{0,1}};
);];


NLSWeqns=Function[{},(
MaxSpeed[h_,hu_]:=Abs[hu/h]+ Sqrt[g h];
Flux[h_,hu_]={hu,hu^2/h+h^2/2};
P[h_,hu_]={0,0}; (*source term from flow vars*)
T[h_,hu_]={0,0}; (*source term from hb*)
dT[h_,hu_]={0,-h g}; (*source term from dx hb*)
A0[h_,hu_]={{0,0},{0,0}};
F[h_,hu_]={{1,0},{0,1}};
);];


DSM0eqns=Function[{},(
MaxSpeed[h_,hu_,hw_,q_]:=Abs[hu/h]+Sqrt[g h];
Flux[h_,hu_,hw_,q_]={hu,hu^2/h+h^2/2+h q,hw hu/h,0};
P[h_,hu_,hw_,q_]={0,0, q, 0}; (*source term from flow vars*)
T[h_,hu_,hw_,q_]={0 ,0,0,0}; (*source term from hb*)
dT[h_,hu_,hw_,q_]={0 ,-q -h g,0,0}; (*source term from dx hb*)
A0[h_,hu_,hw_,q_]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
F[h_,hu_,hw_,q_]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,0}};
);];


DLDSM0eqns=Function[{},(
MaxSpeed[h_,hu_,hw_,q_]:=Abs[hu/h]+1/Fr Sqrt[g h];
Flux[h_,hu_,hw_,q_]={hu,hu^2/h+1/Fr^2 h^2/2+1/Fr^2 h q,hw hu/h,0};
P[h_,hu_,hw_,q_]={0,0, 1/(H^2 Fr^2) q, 0}; (*source term from flow vars*)
T[h_,hu_,hw_,q_]={0 ,0,0,0}; (*source term from hb*)
dT[h_,hu_,hw_,q_]={0 ,-q -h g,0,0}; (*source term from dx hb*)
A0[h_,hu_,hw_,q_]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
F[h_,hu_,hw_,q_]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,0}};
);];


DSM1eqns=Function[{},(
MaxSpeed[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]:=Abs[hu/h]+Sqrt[h+(h\[Alpha]1/h)^2];
Flux[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]={hu,
hu^2/h+1/3h\[Alpha]1^2/h+h^2/2+h q,
2hu h\[Alpha]1/h+h \[Kappa]1,
hu/h hw+ h\[Alpha]1/h h\[Gamma]1 /3,
hw h\[Alpha]1/h+hu/h h\[Gamma]1,
0,
0};
P[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]={0,0,0,(q+\[Kappa]1),3 (-q+\[Kappa]1),0,0}; (*source term from flow vars*)
T[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]={0 ,0,0,0,0,0,0}; (*source term from hb*)
dT[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]={0,-(h+q+\[Kappa]1),0,0,0,0,0}; (*source term from dx hb*)

A0[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]= 
{{0,0,0,0,0,0,0},
{0,0,0,0,0,0,0},
{2 \[Kappa]1,0,hu/h,0,0,0,0},
{0,0,0,0,0,0,0},
{0,0,hw/h,0,0,0,0},
{0,0,0,0,0,0,0},
{0,0,0,0,0,0,0}
};
F[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]=DiagonalMatrix[{1,1,1,1,1,0,0}];
);];


DLDSM1eqns=Function[{},(
MaxSpeed[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]:=Abs[hu/h]+1/Fr Sqrt[g h];
Flux[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]={hu,
hu^2/h+1/3h\[Alpha]1^2/h+1/Fr^2 h^2/2+1/Fr^2 h q,
2hu h\[Alpha]1/h+1/Fr^2 h \[Kappa]1,
hu hw/h+h\[Alpha]1 h\[Gamma]1/(3h),
hw h\[Alpha]1/h+hu h\[Gamma]1/h,
0,
0};
P[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]={0,0,0,1/(H^2 Fr^2)(q+\[Kappa]1),3/(H^2 Fr^2) (-q+\[Kappa]1),0,0}; (*source term from flow vars*)
T[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]={0 ,0,0,0,0,0,0}; (*source term from hb*) (*not scaled...*)
dT[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]={0,-(h+q+\[Kappa]1),0,0,0,0,0}; (*source term from dx hb*)

A0[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]=
{{0,0,0,0,0,0,0},
{0,0,0,0,0,0,0},
{2/Fr^2 \[Kappa]1,0,hu/h,0,0,0,0},
{0,0,0,0,0,0,0},
{0,0,hw/h,0,0,0,0},
{0,0,0,0,0,0,0},
{0,0,0,0,0,0,0}
};
F[h_,hu_,h\[Alpha]1_,hw_,h\[Gamma]1_,q_,\[Kappa]1_]=DiagonalMatrix[{1,1,1,1,1,0,0}];
);];


DLDSM2eqns=Function[{},(
MaxSpeed[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]:=Abs[hu/h]+1/Fr Sqrt[g h] (* max speed probably not correct !! *);
Flux[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]={
hu,
1/Fr^2(h^2/2+h q)+1/(3h)h\[Alpha]1^2+1/(5h)h\[Alpha]2^2+1/h hu^2,
1/Fr^2(h \[Kappa]1)+4/(5h)h\[Alpha]1 h\[Alpha]2+2/h h\[Alpha]1 hu,
1/Fr^2(h \[Kappa]2)+2/(3h)h\[Alpha]1^2+2/(7h)h\[Alpha]2^2+2/h h\[Alpha]2 hu,
1/(3h)h\[Alpha]1 h\[Gamma]1+1/(5h)h\[Alpha]2 h\[Gamma]2 +1/h hu hw,
2/(5h) h\[Alpha]2 h\[Gamma]1+2/(5h) h\[Alpha]1 h\[Gamma]2+1/h h\[Gamma]1 hu +1/h h\[Alpha]1 hw,
2/(3h)h\[Alpha]1 h\[Gamma]1+2/(7h)h\[Alpha]2 h\[Gamma]2+1/h h\[Gamma]2 hu+1/h h\[Alpha]2 hw,
0,
0,
0};
P[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]={0,0,0,0,
1/(H^2 Fr^2) (q+\[Kappa]1+\[Kappa]2),
3/(H^2 Fr^2) (-q+\[Kappa]1+\[Kappa]2),
5/(H^2 Fr^2)(q-\[Kappa]1+\[Kappa]2),
1/h(-hw+h\[Gamma]1-h\[Gamma]2),
3/h(-hw-h\[Gamma]1+h\[Gamma]2),
5/h(-hw-h\[Gamma]1-h\[Gamma]2)}; (*source term from flow vars*)
T[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]=ConstantArray[0,10]; (*source term from hb*)
dT[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]={0,-(h+q+\[Kappa]1+\[Kappa]2),-6\[Kappa]2,0,0,0,0,0,0,0}; (*source term from dx hb*)
A0[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]={
{0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0},
{1/Fr^2 (2\[Kappa]1-3\[Kappa]2),0,-(1/5) h\[Alpha]2/h+hu/h,1/5 h\[Alpha]1/h,0,0,0,0,0,0},
{3/Fr^2 \[Kappa]2,0,h\[Alpha]1/h,1/7 h\[Alpha]2/h+hu/h,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0},
{0,0,-(1/5) h\[Gamma]2/h+hw/h,1/5h\[Gamma]1/h,0,0,0,0,0,0},
{0,0,h\[Gamma]1/h,1/7 h\[Gamma]2/h+hw/h,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0}
};
F[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]=DiagonalMatrix[{1,1,1,1,1,1,1,0,0,0}];
);];


DSM2eqns=Function[{},(
MaxSpeed[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]:=Abs[hu/h]+Sqrt[h+(h\[Alpha]1/h)^2];
Flux[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]={
hu,
h^2/2+h q+1/(3h)h\[Alpha]1^2+1/(5h)h\[Alpha]2^2+1/h hu^2,
h \[Kappa]1+4/(5h)h\[Alpha]1 h\[Alpha]2+2/h h\[Alpha]1 hu,
h \[Kappa]2+2/(3h)h\[Alpha]1^2+2/(7h)h\[Alpha]2^2+2/h h\[Alpha]2 hu,
1/(3h)h\[Alpha]1 h\[Gamma]1+1/(5h)h\[Alpha]2 h\[Gamma]2 +1/h hu hw,
2/(5h) h\[Alpha]2 h\[Gamma]1+2/(5h) h\[Alpha]1 h\[Gamma]2+1/h h\[Gamma]1 hu +1/h h\[Alpha]1 hw,
2/(3h)h\[Alpha]1 h\[Gamma]1+2/(7h)h\[Alpha]2 h\[Gamma]2+1/h h\[Gamma]2 hu+1/h h\[Alpha]2 hw,
0,
0,
0};
P[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]={0,0,0,0,
q+\[Kappa]1+\[Kappa]2,
3 (-q+\[Kappa]1+\[Kappa]2),
5(q-\[Kappa]1+\[Kappa]2),
1/h(-hw+h\[Gamma]1-h\[Gamma]2),
3/h(-hw-h\[Gamma]1+h\[Gamma]2),
5/h(-hw-h\[Gamma]1-h\[Gamma]2)}; (*source term from flow vars*)
T[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]=ConstantArray[0,10]; (*source term from hb*)
dT[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]={0,-(h+q+\[Kappa]1+\[Kappa]2),-6\[Kappa]2,0,0,0,0,0,0,0}; (*source term from dx hb*)
A0[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]={
{0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0},
{2 \[Kappa]1-3\[Kappa]2,0,-(1/5) h\[Alpha]2/h+hu/h,1/5 h\[Alpha]1/h,0,0,0,0,0,0},
{3\[Kappa]2,0,h\[Alpha]1/h,1/7 h\[Alpha]2/h+hu/h,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0},
{0,0,-(1/5) h\[Gamma]2/h+hw/h,1/5h\[Gamma]1/h,0,0,0,0,0,0},
{0,0,h\[Gamma]1/h,1/7 h\[Gamma]2/h+hw/h,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0}
};
F[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,hw_,h\[Gamma]1_,h\[Gamma]2_,q_,\[Kappa]1_,\[Kappa]2_]=DiagonalMatrix[{1,1,1,1,1,1,1,0,0,0}];
);];


DSM3eqns=Function[{},(
MaxSpeed[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,h\[Alpha]3_,hw_,h\[Gamma]1_,h\[Gamma]2_,h\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]:=Abs[hu/h]+Sqrt[ h+(h\[Alpha]1/h)^2];
Flux[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,h\[Alpha]3_,hw_,h\[Gamma]1_,h\[Gamma]2_,h\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]={
hu,
h^2/2+h q+1/(3h)h\[Alpha]1^2+1/(5h)h\[Alpha]2^2+1/(7h)h\[Alpha]3^2+1/h hu^2,
h \[Kappa]1+4/(5h)h\[Alpha]1 h\[Alpha]2+18/(35h)h\[Alpha]2 h\[Alpha]3+2/h h\[Alpha]1 hu,
h \[Kappa]2+2/(3h)h\[Alpha]1^2+2/(7h)h\[Alpha]2^2+2/h h\[Alpha]2 hu+6/(7h)h\[Alpha]1 h\[Alpha]3+4/(21 h)h\[Alpha]3^2,
h \[Kappa]3+6/(5h)h\[Alpha]1 h\[Alpha]2+8/15h h\[Alpha]2 h\[Alpha]3+2/h h\[Alpha]3 hu,
1/(3h)h\[Alpha]1 h\[Gamma]1+1/(5h)h\[Alpha]2 h\[Gamma]2 +1/(7h)h\[Alpha]3 h\[Gamma]3+1/h hu hw,
2/(5h) h\[Alpha]2 h\[Gamma]1+2/(5h) h\[Alpha]1 h\[Gamma]2+9/(35h)h\[Alpha]3 h\[Gamma]2+9/(35h)h\[Alpha]2 h\[Gamma]3+1/h h\[Gamma]1 hu +1/h h\[Alpha]1 hw,
2/(3h)h\[Alpha]1 h\[Gamma]1+3(7h)h\[Alpha]3 h\[Gamma]1+2/(7h)h\[Alpha]2 h\[Gamma]2+3/(7h)h\[Alpha]1 h\[Gamma]3+4/(21h)h\[Alpha]3 h\[Gamma]3+1/h h\[Gamma]2 hu+1/h h\[Alpha]2 hw,
3/(5h)h\[Alpha]2 h\[Gamma]1+3/(5h) h\[Alpha]1 h\[Gamma]2+4/(15h)h\[Alpha]3 h\[Gamma]2+4/15h h\[Alpha]2 h\[Gamma]3+1/h h\[Gamma]3 hu+1/h h\[Alpha]3 hw,
0,
0,
0,
0};
P[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,h\[Alpha]3_,hw_,h\[Gamma]1_,h\[Gamma]2_,h\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]={0,0,0,0,0,
q+\[Kappa]1+\[Kappa]2+\[Kappa]3,
3 (-q+\[Kappa]1+\[Kappa]2+\[Kappa]3),
5(q-\[Kappa]1+\[Kappa]2+\[Kappa]3),
7(-q+\[Kappa]1-\[Kappa]2+\[Kappa]3),
1/h(-hw+h\[Gamma]1-h\[Gamma]2+h\[Gamma]3),
3/h(-hw-h\[Gamma]1+h\[Gamma]2-h\[Gamma]3),
5/h(-hw-h\[Gamma]1-h\[Gamma]2+h\[Gamma]3),
7/h(-hw-h\[Gamma]1-h\[Gamma]2-h\[Gamma]3)}; (*source term from flow vars*)
T[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,h\[Alpha]3_,hw_,h\[Gamma]1_,h\[Gamma]2_,h\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]=ConstantArray[0,13]; (*source term from hb*)
dT[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,h\[Alpha]3_,hw_,h\[Gamma]1_,h\[Gamma]2_,h\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]={0,-(h+q+\[Kappa]1+\[Kappa]2+\[Kappa]3),-6\[Kappa]2,-10\[Kappa]3,0,0,0,0,0,0,0,0,0}; (*source term from dx hb*)
A0[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,h\[Alpha]3_,hw_,h\[Gamma]1_,h\[Gamma]2_,h\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]={
{0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0},
{2 \[Kappa]1-3\[Kappa]2+3\[Kappa]3,0,-(1/5) h\[Alpha]2/h+hu/h,1/5 h\[Alpha]1/h-3/35 h\[Alpha]3/h,3/35 h\[Alpha]2/h,0,0,0,0,0,0,0,0},
{-5 \[Kappa]3+3\[Kappa]2,0,h\[Alpha]1/h-3/7h\[Alpha]3/h,1/7 h\[Alpha]2/h+hu/h,2/7h\[Alpha]1/h+1/21 h\[Alpha]3/h,0,0,0,0,0,0,0,0},
{4\[Kappa]3,0,6/5h\[Alpha]2/h,4/5 h\[Alpha]1/h+2/15 h\[Alpha]3/h,1/5 h\[Alpha]2/h+hu/h,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,-(1/5) h\[Gamma]2/h+hw/h,1/5h\[Gamma]1/h-3/35 h\[Gamma]3/h,3/35h\[Gamma]2/h,0,0,0,0,0,0,0,0},
{0,0,h\[Gamma]1/h-3/7 h\[Gamma]3/h,1/7 h\[Gamma]2/h+hw/h,2/7h\[Gamma]1/h+1/21 h\[Gamma]3/h,0,0,0,0,0,0,0,0},
{0,0,6/5 h\[Gamma]2/h,4/5 h\[Gamma]1/h+2/15 h\[Gamma]3/h,1/5 h\[Gamma]2/h+hw/h,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0}
};
F[h_,hu_,h\[Alpha]1_,h\[Alpha]2_,h\[Alpha]3_,hw_,h\[Gamma]1_,h\[Gamma]2_,h\[Gamma]3_,q_,\[Kappa]1_,\[Kappa]2_,\[Kappa]3_]=DiagonalMatrix[{1,1,1,1,1,1,1,1,1,0,0,0,0}];
);];
