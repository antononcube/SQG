(* ::Package:: *)

(* ::Input:: *)
(* ==========================SQG2.nb==========================Lorentzian (Minkowski),covariant equal-trace stepper.-mixed Lorentzian't Hooft (\[Mu]\[DownArrow],\[Nu]\[UpArrow])-Urbantke g(P) used for all spacetime contractions-S(P)=P g(P) P^T (internal),\[CapitalOmega] mixed (spacetime),H internal-Duality: \bar\[Eta]^i_{\[Mu]}{}^{\[Nu]} tr(\[Sigma]_i P^{\[Mu]} \[Delta]P_{\[Nu]})-Traces:tr \[CapitalOmega]=tr H=u (metric-free)============================================================*)

(*tolerances/precision*)
If[!ValueQ[$SQGTol],$SQGTol=1.*^-12];
If[!ValueQ[$SQGWP],$SQGWP=80];

(*-----Pauli& internal basis T_A-----*)
pauli={{{0,1},{1,0}},{{0,-I},{I,0}},{{1,0},{0,-1}}};
tBasis[A_]:=If[A==1,I IdentityMatrix[2],pauli[[A-1]]];

(*-----Lorentzian mixed't Hooft (\[Mu]\[DownArrow],\[Nu]\[UpArrow]),\[CurlyEpsilon]^{0123}=+1;slots 1..4 map 0..3-----*)
Clear[EtaMixed,BarEtaMixed];
EtaMixed:=EtaMixed=Module[{e3=LeviCivitaTensor[3],\[Delta]=KroneckerDelta,E},E=ConstantArray[0.,{3,4,4}];
Do[E[[i,j+1,k+1]]=e3[[i,j,k]],{i,1,3},{j,1,3},{k,1,3}];
Do[E[[i,j+1,1]]=+I \[Delta][i,j],{i,1,3},{j,1,3}];
Do[E[[i,1,j+1]]=-I \[Delta][i,j],{i,1,3},{j,1,3}];
Developer`ToPackedArray@E];
BarEtaMixed:=BarEtaMixed=Module[{e3=LeviCivitaTensor[3],\[Delta]=KroneckerDelta,E},E=ConstantArray[0.,{3,4,4}];
Do[E[[i,j+1,k+1]]=e3[[i,j,k]],{i,1,3},{j,1,3},{k,1,3}];
Do[E[[i,j+1,1]]=-I \[Delta][i,j],{i,1,3},{j,1,3}];
Do[E[[i,1,j+1]]=+I \[Delta][i,j],{i,1,3},{j,1,3}];
Developer`ToPackedArray@E];

(*-----spacetime helpers-----*)
Pmu2x2[P_,\[Mu]_Integer?Positive]:=Sum[tBasis[A]*P[[A,\[Mu]]],{A,1,4}];

(*LOWER only (P stored with upper \[Mu]):P^A{}_\[Nu]=g_{\[Nu]\[Mu]} P^{A\[Mu]}*)
Clear[LowerP,gDot];
LowerP[P_,gDn_]:=P . Transpose[gDn];
gDot[nDn_,vUp_,gDn_]:=nDn . gDn . vUp;  (*in case you need n^\[Mu] g_{\[Mu]\[Nu]} v^\[Nu]*)

(*-----curvature F from P (for Urbantke)-----*)
Fhat[P_]:=Module[{F=ConstantArray[0.,{3,4,4}],Pm,Pn,comm},Do[Pm=Pmu2x2[P,\[Mu]];Pn=Pmu2x2[P,\[Nu]];comm=Pm . Pn-Pn . Pm;
Do[F[[i,\[Mu],\[Nu]]]=Re[-I*Tr[pauli[[i]] . comm]];
F[[i,\[Nu],\[Mu]]]=-F[[i,\[Mu],\[Nu]]],{i,1,3}],{\[Mu],1,4},{\[Nu],\[Mu]+1,4}];
Developer`ToPackedArray[F]];

(*-----Urbantke g(P),det-normalized-----*)
Clear[UrbantkeMetric];
UrbantkeMetric[P_,tol_:$SQGTol]:=Module[{F = FfromP[],FF, G},
(*	FF ={F[[2]] . F[[3]]- F[[3]] . F[[2]],F[[3]] . F[[1]]- F[[1]] . F[[3]],F[[1]] . F[[2]]- F[[2]] . F[[1]]};  G =-1/6  Sum[ FF[[k]] . F[[k]] ,{k,3}];*)
detG=Det@N@G;
gDn=If[NumericQ[detG]&&detG=!=0,N@Chop[G/Abs[detG]^(1/4),tol],IdentityMatrix[4]];
gUp=Inverse[gDn];
<|"g"->gDn,"ginv"->gUp|>];

(*-----covariant internal Gram and A(\[CapitalOmega])-----*)
Clear[Sgram,AfromOmega,Skew];
Sgram[P_,g_,tol_:$SQGTol]:=Module[{},Developer`ToPackedArray@N@Chop[P . g . Transpose[P],tol]];
AfromOmega[P_,g_,Om_,tol_:$SQGTol]:=Module[{},Developer`ToPackedArray@N@Chop[P . g . Om . Transpose[P],tol]];
Skew[M_]:=(M-Transpose[M])/2.;

(*-----stable internal 4\[Times]4 solve for Xi-----*)
Clear[XiFrom];
XiFrom[P_,g_,S_,Om_,H_,mu_:Automatic]:=Module[{A=AfromOmega[P,g,Om],\[Mu],wp=$SQGWP,SS,RHS,X},\[Mu]=If[mu===Automatic,1.*^-8,mu];
SS=SetPrecision[N@Chop[S,$SQGTol],wp]+SetPrecision[\[Mu],wp] IdentityMatrix[4];
RHS=SetPrecision[N@Chop[Skew[A],$SQGTol],wp];
X=LinearSolve[SS,RHS];
-X+H];

(*-----duality residuals: \bar\[Eta]^i_{\[Mu]}{}^{\[Nu]} tr(\[Sigma]_i P^{\[Mu]} \[Delta]P_{\[Nu]})-----*)
Clear[DualityBlocksFor];
DualityBlocksFor[P_,g_,dP_,tol_:$SQGTol]:=Module[
{be=BarEtaMixed,r1=ConstantArray[0.,{3,3}],r2=ConstantArray[0.,{3,3}],Be, tAB, PdP, PdM,dPP,MdP},
(* Be = {g . be[[1]],g . be[[2]],g . be[[3]]};*)
(* (3,4,4) *)
tAB = Table[tBasis[A] . tBasis[B],{A,4},{B,4}];
(* (4,4,2,2)*)
PdP = Table[P . Be[[a]] . Transpose[dP],{a,3}];
(* (3,4,4)*)
dPP = Table[dP . Be[[a]] . Transpose[P],{a,3}];
(* (3,4,4)*)
(* MdP =TensorContract[TensorProduct[ PdP,tAB],{{2,1},{3,2}}];*)
(* (3,2,2) *)
PdM =TensorContract[TensorProduct[ dPP,tAB],{{2,1},{3,2}}];
(* (3,2,2) *)
Do[
r1[[a,i]]+=Im@Tr[pauli[[i]] . PdM[[a]]];
(*         r2[[a,i]]+=Im@Tr[pauli[[i]] . MdP[[a]]],*)
{a,1,3},{i,1,3}];
Developer`ToPackedArray@Join[Flatten[r1],Flatten[r2]]];

(*-----internal projectors \[CapitalPi]_k (trace=1) from S-----*)
Clear[InternalProjectors];
InternalProjectors[S_,tol_:$SQGTol]:=Module[{vals,vecs,V},{vals,vecs}=Eigensystem[S];
V=Transpose@Orthogonalize[vecs];
V=N@Chop[V,tol];
Table[With[{v=Normalize[V[[All,k]]]},Developer`ToPackedArray@Outer[v,v]],{k,1,4}]];

(*-----build 20\[Times]20 (18 duality+metric-free tr\[CapitalOmega],trH)-----*)
Clear[BuildSystem20];
BuildSystem20[P_,g_,S_,u_?NumericQ,tol_:$SQGTol]:=Module[{Pnum=N@Chop[P,tol],\[CapitalPi]s,cols,Oe,He,dP,L18,trO,trH,L,b},
\[CapitalPi]s=InternalProjectors[S,tol];
cols=ConstantArray[0.,{20,18}];
Do[Oe=ConstantArray[0.,{4,4}];Oe[[ai,bj]]=1.;
He=ConstantArray[0.,{4,4}];
dP=Oe . Pnum-Pnum . XiFrom[Pnum,g,S,Oe,He,Automatic];
cols[[(ai-1)*4+bj,;;]]=DualityBlocksFor[Pnum,g,dP,tol],{ai,1,4},{bj,1,4}];
Do[Oe=ConstantArray[0.,{4,4}];
He=\[CapitalPi]s[[k]];
dP=-Pnum . He;
cols[[16+k,;;]]=DualityBlocksFor[Pnum,g,dP,tol],{k,1,4}];
L18=Transpose[cols];
trO=Join[Flatten@IdentityMatrix[4],ConstantArray[0.,4]];
trH=Join[ConstantArray[0.,16],ConstantArray[1.,4]];
L=Join[L18,{trO,trH}];
b=Join[ConstantArray[0.,18],{u,u}];
{Developer`ToPackedArray[L],Developer`ToPackedArray[b],\[CapitalPi]s}];

(*-----optional one-shot step consistency checker-----*)
Clear[AssertStepConsistency];
AssertStepConsistency[P_,Om_,H_,tol_:$SQGTol]:=Module[{ok=True,msg={}},If[Dimensions[P]=!={4,4},ok=False;AppendTo[msg,"P not 4x4"]];
If[Dimensions[Om]=!={4,4},ok=False;AppendTo[msg,"\[CapitalOmega] not 4x4 (spacetime)"]];
If[Dimensions[H]=!={4,4},ok=False;AppendTo[msg,"H not 4x4 (internal)"]];
If[Abs[Tr[Om]-Tr[H]]>10 tol,ok=False;AppendTo[msg,"tr\[CapitalOmega] != trH: "<>ToString@NumberForm[{Tr[Om],Tr[H]},{8,3}]]];
If[!ok,Print["[step-check] ",StringRiffle[msg,"; "]]];
ok];

(*-----solve 20\[Times]20,single SDE step;enforce P^0 imaginary;assert once-----*)
If[!ValueQ[$AssertStepOnce],$AssertStepOnce=True];

Clear[SolveOmegaH20];
SolveOmegaH20[P_,g_,S0_,u_,tol_:$SQGTol,wp_:$SQGWP]:=Module[{L,b,\[CapitalPi]s,U,S,V,s,smax,Sinv,v},{L,b,\[CapitalPi]s}=BuildSystem20[P,g,S0,u,tol];
{U,S,V}=SingularValueDecomposition[SetPrecision[L,wp]];
s=Diagonal[S];smax=Max[s,0.];
Sinv=DiagonalMatrix@Map[If[#>tol*smax,1./#,0.]&,s];
v=V . Sinv . Transpose[U] . SetPrecision[b,wp];
<|"\[CapitalOmega]"->ArrayReshape[v[[1;;16]],{4,4}],"h"->v[[17;;20]],"\[CapitalPi]s"->\[CapitalPi]s|>];

Clear[StepSDE20];
StepSDE20[P_,u_,d\[Theta]_:1.,stabEvery_:8,stepIndex_:1,tol_:$SQGTol,wp_:$SQGWP]:=Module[{g= UrbantkeMetric[P,tol],S,sol,Om,h,\[CapitalPi]s,Hint,Xi,dP,Pnew,det},
S=Sgram[P,g,tol];
sol=SolveOmegaH20[P,g,S,u,tol,wp];
Om=sol["\[CapitalOmega]"];h=sol["h"];\[CapitalPi]s=sol["\[CapitalPi]s"];
Hint=Sum[h[[k]] \[CapitalPi]s[[k]],{k,1,4}];(*internal H*)If[$AssertStepOnce&&stepIndex==1,AssertStepConsistency[P,Om,Hint,$SQGTol];$AssertStepOnce=False];
Xi=XiFrom[P,g,S,Om,Hint,Automatic];(*internal 4\[Times]4*)dP=Om . P-P . Xi;(*mixed spacetime/internal*)Pnew=P+d\[Theta] dP;
If[Mod[stepIndex,stabEvery]==0,det=Det[Pnew];
If[det<0,Pnew[[1,All]]*=-1;det=Det[Pnew]];
Pnew=Pnew/Abs[det]^(1/4)];
(*enforce abelian row purely imaginary*)Pnew[[1,All]]=I*Re[-I*Pnew[[1,All]]];
(*optional:keep non-abelian rows real*)(*Pnew[[2;;4,All]]=Re[Pnew[[2;;4,All]]];*)<|"P"->Developer`ToPackedArray[Pnew],"\[CapitalOmega]"->Om,"h"->h|>];

(*-----invariants (same API)-----*)
SafeSqrtSym[X_,tol_:$SQGTol]:=Module[{S=Symmetrize[X],vals,vecs,vclip},{vals,vecs}=Eigensystem[S];
vclip=vals/. x_/;x<0&&Abs[x]<tol->0.;
vecs . DiagonalMatrix[Sqrt[Clip[vclip,{0,\[Infinity]}]]] . Transpose[vecs]];

XfromP[P_,g_]:=Module[{F=Fhat[P],X=ConstantArray[0.,{3,3}], TF},
TF = {g . F[[1]] . g,g . F[[2]] . g,g . F[[3]] . g};
Do[
X[[j,i]]=X[[i,j]]= 1/2TF[[i]] . F[[j]],
{i,1,3},{j,i,3}];
Developer`ToPackedArray[X]];

CosmoInvariants[P_,g_,z_?NumericQ,s_:+1]:=Module[{X,Y,trY,R,W},X=XfromP[P,g];Y=SafeSqrtSym[X];trY=Tr[Y];
R=4.*s*trY;W=z*R;
<|"trY"->trY,"R"->R,"W"->W|>];

WfromP[P_,z_,s_:+1]:=CosmoInvariants[P,z,s]["W"];




(* ====================== Euclidean patch (appended) ======================
   Patched on 2025-09-18T05:25:29.686775 by Eily
   - Switch to Euclidean tetrad basis throughout SQG2 (no Lorentzian *i factors).
   - EtaMixed/BarEtaMixed become real Euclidean 't Hooft symbols.
   - UrbantkeMetric uses the fast triple-product (for optional checks only).
   - Default metric g is δ (Identity); SDE code can keep passing g unchanged.
========================================================================= *)

(* Euclidean 't Hooft symbols η^i_{AB} (self-dual) and \barη (anti-self-dual).
   Indices A,B=1..4 map (r̂, θ̂, φ̂, t_E) if you prefer; here we keep slots 1..4.
   Conventions: η^i_{ab} = ε_{iab},  η^i_{a4} = +δ_{ia},  η^i_{4a} = -δ_{ia}. *)
Clear[EtaMixed, BarEtaMixed];
EtaMixed := EtaMixed = Module[{e3 = LeviCivitaTensor[3], δ = KroneckerDelta, E},
  E = ConstantArray[0., {3, 4, 4}];
  Do[E[[i, j, k]] = e3[[i, j, k]], {i, 1, 3}, {j, 1, 3}, {k, 1, 3}];
  Do[E[[i, j, 4]] = +δ[i, j], {i, 1, 3}, {j, 1, 3}];
  Do[E[[i, 4, j]] = -δ[i, j], {i, 1, 3}, {j, 1, 3}];
  Developer`ToPackedArray @ E
];
BarEtaMixed := BarEtaMixed = Module[{e3 = LeviCivitaTensor[3], δ = KroneckerDelta, E},
  E = ConstantArray[0., {3, 4, 4}];
  Do[E[[i, j, k]] = e3[[i, j, k]], {i, 1, 3}, {j, 1, 3}, {k, 1, 3}];
  Do[E[[i, j, 4]] = -δ[i, j], {i, 1, 3}, {j, 1, 3}];
  Do[E[[i, 4, j]] = +δ[i, j], {i, 1, 3}, {j, 1, 3}];
  Developer`ToPackedArray @ E
];

(* P→2×2 kept as-is; treat columns as tetrad indices. *)

(* Fast Euclidean Urbantke (optional; not used in SDE) *)
Clear[UrbantkeMetric];
UrbantkeMetric[P_, tol_ : $SQGTol] := Module[{F, FF1, FF2, FF3, G, detG, gDn, gUp},
  F = Fhat[P];
  FF1 = F[[2]].F[[3]] - F[[3]].F[[2]];
  FF2 = F[[3]].F[[1]] - F[[1]].F[[3]];
  FF3 = F[[1]].F[[2]] - F[[2]].F[[1]];
  G = -(FF1.F[[1]] + FF2.F[[2]] + FF3.F[[3]])/6.;
  G = Chop[(G + Transpose[G])/2, tol];
  detG = Det @ N @ G;
  gDn = If[NumericQ[detG] && detG =!= 0, N @ Chop[G/Abs[detG]^(1/4), tol], IdentityMatrix[4]];
  gUp = Inverse[gDn];
  <|"g" -> gDn, "ginv" -> gUp, "densitized" -> G|>
];

(* Euclidean metric helper *)
If[!ValueQ[$EuclidMetric], $EuclidMetric = IdentityMatrix[4]];
