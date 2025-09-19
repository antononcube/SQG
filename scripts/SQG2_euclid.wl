(* ::Package:: *)

(* ::Input:: *)
(* ========================== SQG2.wl (Euclidean tetrad, metric-free) ========================== *)

(*tolerances/precision*)
If[!ValueQ[$SQGTol],$SQGTol=1.*^-12];
If[!ValueQ[$SQGWP],$SQGWP=80];

(*-----Pauli& internal basis T_A-----*)
pauli={{{0,1},{1,0}},{{0,-I},{I,0}},{{1,0},{0,-1}}};
tBasis[A_]:=If[A==1,I IdentityMatrix[2],pauli[[A-1]]];



(* ----- Euclidean mixed 't Hooft symbols (time slot = 1) ----- *)
(* A,B = 1..4 with 1 ≡ t_E; spatial slots 2..4 map to a=1..3.
   η^i_{ab} = ε_{iab},   η^i_{a,1} = +δ_{ia},   η^i_{1,a} = -δ_{ia}
   arη^i_{ab} = ε_{iab}, arη^i_{a,1} = -δ_{ia}, arη^i_{1,a} = +δ_{ia} *)
Clear[EtaMixed, BarEtaMixed];
EtaMixed := EtaMixed = Module[{e3 = LeviCivitaTensor[3], δ = KroneckerDelta, E},
  E = ConstantArray[0., {3, 4, 4}];
  Do[E[[i, a+1, b+1]] = e3[[i, a, b]], {i, 1, 3}, {a, 1, 3}, {b, 1, 3}];
  Do[E[[i, a+1, 1]] =  δ[i, a],       {i, 1, 3}, {a, 1, 3}];
  Do[E[[i, 1, a+1]] = -δ[i, a],       {i, 1, 3}, {a, 1, 3}];
  Developer`ToPackedArray @ E
];
BarEtaMixed := BarEtaMixed = Module[{e3 = LeviCivitaTensor[3], δ = KroneckerDelta, E},
  E = ConstantArray[0., {3, 4, 4}];
  Do[E[[i, a+1, b+1]] = e3[[i, a, b]], {i, 1, 3}, {a, 1, 3}, {b, 1, 3}];
  Do[E[[i, a+1, 1]] = -δ[i, a],       {i, 1, 3}, {a, 1, 3}];
  Do[E[[i, 1, a+1]] =  δ[i, a],       {i, 1, 3}, {a, 1, 3}];
  Developer`ToPackedArray @ E
];


(*-----spacetime helpers-----*)
Pmu2x2[P_,\[Mu]_Integer?Positive]:=Sum[tBasis[A]*P[[A,\[Mu]]],{A,1,4}];

(*LOWER only (P stored with upper \[Mu]):P^A{}_\[Nu]=g_{\[Nu]\[Mu]} P^{A\[Mu]}*)

(*-----curvature F from P (for Urbantke)-----*)
Fhat[P_]:=Module[{F=ConstantArray[0.,{3,4,4}],Pm,Pn,comm},Do[Pm=Pmu2x2[P,\[Mu]];Pn=Pmu2x2[P,\[Nu]];comm=Pm . Pn-Pn . Pm;
Do[F[[i,\[Mu],\[Nu]]]=Re[-I*Tr[pauli[[i]] . comm]];
F[[i,\[Nu],\[Mu]]]=-F[[i,\[Mu],\[Nu]]],{i,1,3}],{\[Mu],1,4},{\[Nu],\[Mu]+1,4}];
Developer`ToPackedArray[F]];

(* Clear[UrbantkeMetric];
UrbantkeMetric[P_, tol_:$SQGTol] := Module[{F, FF1, FF2, FF3, G, d, gDn},
  F = Fhat[P];                                (* F is 3×4×4 of antisymmetric 2-forms *)
  FF1 = F[[2]].F[[3]] - F[[3]].F[[2]];        (* ε_{23k} F^2 F^3, etc. *)
  FF2 = F[[3]].F[[1]] - F[[1]].F[[3]];
  FF3 = F[[1]].F[[2]] - F[[2]].F[[1]];
  G   = -(FF1.F[[1]] + FF2.F[[2]] + FF3.F[[3]])/6.;   (* √g g = -(1/6) ε F F F *)
  G   = Chop[(G + Transpose[G])/2, tol];              (* symmetric numerically *)
  d   = Det @ N @ G;
  gDn = If[NumericQ[d] && d =!= 0, N @ Chop[G/Abs[d]^(1/4), tol], IdentityMatrix[4]];
  <|"densitized" -> G, "g" -> gDn, "ginv" -> Inverse[gDn]|>
]; *)

(*-----covariant internal Gram and A(\[CapitalOmega])-----*)
Clear[Sgram,AfromOmega,Skew];
Sgram[P_,tol_:$SQGTol]:=Module[{},Developer`ToPackedArray@N@Chop[P . Transpose[P],tol]];
AfromOmega[P_,Om_,tol_:$SQGTol]:=Module[{},Developer`ToPackedArray@N@Chop[P . Om . Transpose[P],tol]];
Skew[M_]:=(M-Transpose[M])/2.;

(*-----stable internal 4\[Times]4 solve for Xi-----*)
Clear[XiFrom];
XiFrom[P_,S_,Om_,H_,mu_:Automatic]:=Module[{A=AfromOmega[P,Om],\[Mu],wp=$SQGWP,SS,RHS,X},\[Mu]=If[mu===Automatic,1.*^-8,mu];
SS=SetPrecision[N@Chop[S,$SQGTol],wp]+SetPrecision[\[Mu],wp] IdentityMatrix[4];
RHS=SetPrecision[N@Chop[Skew[A],$SQGTol],wp];
X=LinearSolve[SS,RHS];
-X+H];

(*-----duality residuals: \bar\[Eta]^i_{\[Mu]}{}^{\[Nu]} tr(\[Sigma]_i P^{\[Mu]} \[Delta]P_{\[Nu]})-----*)
Clear[DualityBlocksFor];
DualityBlocksFor[P_,dP_,tol_:$SQGTol]:=Module[{
    be=BarEtaMixed,
    r1=ConstantArray[0.,{3,3}],
    r2=ConstantArray[0.,{3,3}],
    tAB, PdP, PdM, dPP, MdP
  },
  tAB = Table[tBasis[A] . tBasis[B],{A,4},{B,4}];
  (* (4,4,2,2)*)
  PdP = Table[P . be[[a]] . Transpose[dP],{a,3}];
  (* (3,4,4)*)
  dPP = Table[dP . be[[a]] . Transpose[P],{a,3}];
  (* (3,2,2) *)                          (* (3,4,4)(4,4,2,2)*)
  MdP =TensorContract[TensorProduct[ PdP,tAB],{{2,4},{3,5}}];
  (* (3,2,2) *)                          (* (3,4,4)(4,4,2,2)*)
  PdM =TensorContract[TensorProduct[ dPP,tAB],{{2,4},{3,5}}];
  (* (3,2,2) *)
  Do[
    r1[[a,i]]+=Im@Tr[pauli[[i]] . PdM[[a]]];
    r2[[a,i]]+=Im@Tr[pauli[[i]] . MdP[[a]]],
    {a,1,3}, {i,1,3}
  ];
  Developer`ToPackedArray@Join[Flatten[r1],Flatten[r2]]
];

(*-----internal projectors \[CapitalPi]_k (trace=1) from S-----*)
Clear[InternalProjectors];
InternalProjectors[S_, tol_:$SQGTol] := Module[{vals, vecs, V},
  {vals,vecs}=Eigensystem[S];
  Table[
    With[{v = vecs[[k]]},
      Developer`ToPackedArray@TensorProduct[v,v]
    ],
    {k,1,4}
  ]
];

(*-----build 20\[Times]20 (18 duality+metric-free tr\[CapitalOmega],trH)-----*)
Clear[BuildSystem20];
BuildSystem20[P_,S_,u_?NumericQ,tol_:$SQGTol]:=Module[{Pnum=N@Chop[P,tol],\[CapitalPi]s,cols,Oe,He,dP,L18,trO,trH,L,b},
\[CapitalPi]s=InternalProjectors[S,tol];
cols=ConstantArray[0.,{20,18}];
Do[Oe=ConstantArray[0.,{4,4}];Oe[[ai,bj]]=1.;
He=ConstantArray[0.,{4,4}];
dP=Oe . Pnum-Pnum . XiFrom[Pnum,S,Oe,He,Automatic];
cols[[(ai-1)*4+bj,;;]]=DualityBlocksFor[Pnum,dP,tol],{ai,1,4},{bj,1,4}];
Do[Oe=ConstantArray[0.,{4,4}];
He=\[CapitalPi]s[[k]];
dP=-Pnum . He;
cols[[16+k,;;]]=DualityBlocksFor[Pnum,dP,tol],{k,1,4}];
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
SolveOmegaH20[P_, S0_, u_, f_, tol_ : $SQGTol, wp_ : $SQGWP] :=
 Module[{L, b, Πs, U, S, V, s, smax, d, Sinv, v, i0, zIdx, fUse},
  {L, b, Πs} = BuildSystem20[P, S0, u, tol];

  {U, S, V} = SingularValueDecomposition[SetPrecision[L, wp]];
  s    = Diagonal[S];
  smax = Max[s, 0.];
  Print["after SVD: s=", N[s]];

  (* pseudoinverse diag: nonzeros -> 1/s, zeros -> 0 *)
  d    = If[# > tol*smax, 1./#, 0.] & /@ s;
  Sinv = DiagonalMatrix[d];
  v    = V . Sinv . Transpose[U] . SetPrecision[b, wp];   (* least-norm *)

  (* indices of zero singulars: tail starting at first zero of d *)
  i0   = FirstPosition[d, 0., Missing["NotFound"]];
  zIdx = If[MissingQ[i0], {}, Range[i0[[1]], Length@d]];

  (* match kick length to nullspace dim, or skip if none *)
  fUse = If[zIdx === {}, {}, Take[PadRight[Flatten@{f}, Length@zIdx, 0.], Length@zIdx]];
  If[fUse =!= {}, v = v + V[[All, zIdx]] . fUse];

  <|
    "\[CapitalOmega]" -> ArrayReshape[v[[1 ;; 16]], {4, 4}],
    "h"               -> v[[17 ;; 20]],
    "\[CapitalPi]s"   -> Πs
  |>
]


Clear[StepSDE20];
StepSDE20[P_,u_,f_,d\[Theta]_:1.,stabEvery_:8,stepIndex_:1,tol_:$SQGTol,wp_:$SQGWP]:=
Module[{S,sol,Om,h,\[CapitalPi]s,Hint,Xi,dP,Pnew,det},
  S=Sgram[P,tol];
  sol=SolveOmegaH20[P,S,u,f,tol,wp];
  Om=sol["\[CapitalOmega]"];h=sol["h"];\[CapitalPi]s=sol["\[CapitalPi]s"];
  Hint=Sum[h[[k]] \[CapitalPi]s[[k]],{k,1,4}];(*internal H*)If[$AssertStepOnce&&stepIndex==1,AssertStepConsistency[P,Om,Hint,$SQGTol];$AssertStepOnce=False];
  Xi=XiFrom[P,S,Om,Hint,Automatic];(*internal 4\[Times]4*)dP=Om . P-P . Xi;(*mixed spacetime/internal*)Pnew=P+d\[Theta] dP;
  det=Det[Pnew];
  Assert[Im[det]>0 && Abs[Re[det]] < 0.1 tol];
  Pnew=Pnew/Abs[det]^(1/4);
  (*enforce abelian row purely imaginary*)
  Pnew[[1,All]]=I Im[Pnew[[1,All]]];
  (*optional:keep non-abelian rows real*)
  (*Pnew[[2;;4,All]]=Re[Pnew[[2;;4,All]]];*)
  <|"P"->Developer`ToPackedArray[Pnew],"\[CapitalOmega]"->Om,"h"->h|>
];

(*-----invariants (same API)-----*)
SafeSqrtSym[X_,tol_:$SQGTol]:=
Module[{S=Symmetrize[X],vals,vecs,vclip},
{vals,vecs}=Eigensystem[S];
Print["x_vals=",vals];
vclip=vals/. x_/;x<0&&Abs[x]<tol->0.;
vecs . DiagonalMatrix[Sqrt[Clip[vclip,{0,\[Infinity]}]]] . Transpose[vecs]];

XfromP[P_]:=Module[{F=Fhat[P],X=ConstantArray[0.,{3,3}]},
Do[
X[[j,i]]=X[[i,j]]= -1/2Tr[F[[i]] . F[[j]]],
{i,1,3},{j,i,3}];
Developer`ToPackedArray[X]];

CosmoInvariants[P_,z_?NumericQ,s_:+1]:=
Module[{X,Y,trY,R,W},X=XfromP[P];
Y=SafeSqrtSym[X];
trY=Tr[Y];
R=4.*s*trY;W=z*R;
<|"trY"->trY,"R"->R,"W"->W|>];

RfromP[P_,z_,s_:+1]:=CosmoInvariants[P,z,s]["R"];
WfromP[P_,z_,s_:+1]:=CosmoInvariants[P,z,s]["W"];
