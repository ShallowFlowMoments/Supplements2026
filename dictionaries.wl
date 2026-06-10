(* ::Package:: *)

dictLDSM0[Uvec_]:={
u[i_]->Indexed[Uvec[[2]],i],
w[i_]->Indexed[Uvec[[3]],i],
ltag->TRUE}


dictLDSM1[Uvec_]:={
u[i_]->Indexed[Uvec[[2]],i],
a1[i_]->Indexed[Uvec[[3]],i],
w[i_]->Indexed[Uvec[[4]],i],
g1[i_]->Indexed[Uvec[[5]],i],
ltag->True}


dictLDSM2[Uvec_]:={u[i_]->Indexed[Uvec[[2]],i],
a1[i_]->Indexed[Uvec[[3]],i],
a2[i_]->Indexed[Uvec[[4]],i],
w[i_]->Indexed[Uvec[[5]],i],
g1[i_]->Indexed[Uvec[[6]],i],
g2[i_]->Indexed[Uvec[[7]],i],
ltag->True}


dictLDSM3[Uvec_]:={u[i_]->Indexed[Uvec[[2]],i],
a1[i_]->Indexed[Uvec[[3]],i],
a2[i_]->Indexed[Uvec[[4]],i],
a3[i_]->Indexed[Uvec[[5]],i],
w[i_]->Indexed[Uvec[[6]],i],
g1[i_]->Indexed[Uvec[[7]],i],
g2[i_]->Indexed[Uvec[[8]],i],
g3[i_]->Indexed[Uvec[[9]],i],
ltag->True}


dictDSM0[Uvec_]:={
u[i_]->Indexed[Uvec[[2]]/Uvec[[1]],i],
w[i_]->Indexed[Uvec[[3]]/Uvec[[1]],i],
ltag->False}


dictDSM1[Uvec_]:={
u[i_]->Indexed[Uvec[[2]]/Uvec[[1]],i],
a1[i_]->Indexed[Uvec[[3]]/Uvec[[1]],i],
w[i_]->Indexed[Uvec[[4]]/Uvec[[1]],i],
g1[i_]->Indexed[Uvec[[5]]/Uvec[[1]],i],
ltag->False}


dictDSM2[Uvec_]:={u[i_]->Indexed[Uvec[[2]]/Uvec[[1]],i],
a1[i_]->Indexed[Uvec[[3]]/Uvec[[1]],i],
a2[i_]->Indexed[Uvec[[4]]/Uvec[[1]],i],
w[i_]->Indexed[Uvec[[5]]/Uvec[[1]],i],
g1[i_]->Indexed[Uvec[[6]]/Uvec[[1]],i],
g2[i_]->Indexed[Uvec[[7]]/Uvec[[1]],i],
ltag->False}


dictDSM3[Uvec_]:={u[i_]->Indexed[Uvec[[2]]/Uvec[[1]],i],
a1[i_]->Indexed[Uvec[[3]]/Uvec[[1]],i],
a2[i_]->Indexed[Uvec[[4]]/Uvec[[1]],i],
a3[i_]->Indexed[Uvec[[5]]/Uvec[[1]],i],
w[i_]->Indexed[Uvec[[6]]/Uvec[[1]],i],
g1[i_]->Indexed[Uvec[[7]]/Uvec[[1]],i],
g2[i_]->Indexed[Uvec[[8]]/Uvec[[1]],i],
g3[i_]->Indexed[Uvec[[9]]/Uvec[[1]],i],
ltag->False}
