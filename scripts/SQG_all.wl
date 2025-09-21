(* ::Package:: *)

$HistoryLength = 0;
Clear[SQGPrint];
$SQGDiagnostics = True;
If[!ValueQ[$SQGDiagnostics], $SQGDiagnostics = False];
SQGPrint[args___] := If[TrueQ[$SQGDiagnostics], Print[args]];


(* ::Package:: *)
(**)


(* ====================== SQG1.wl (Euclidean tetrad seed) ======================
   Patched on 2025-09-18T05:25:29.686775 by Eily: switch to Euclidean tetrad basis.
   - Work entirely in tetrad indices M=1..4 \[Congruent] (t_E, r\:0302, \[Theta]\:0302, \[CurlyPhi]\:0302), with \[Delta]_{AB}.
   - SeedFromZ now returns SU(2) seed P (4\[Times]4) in tetrad basis; Abelian row left 0.
   - Optional Urbantke check included (Euclidean triple-product form), local-only.
============================================================================= *)

If[!ValueQ[$SQGTol], $SQGTol = 1.*^-12];

(* Curvature scale \[Psi](z) with z = r^3 / r0 (M=1) *)
Clear[psiFromZ];
psiFromZ[z_?NumericQ] := 1./(2. z);

(* Euclidean BH seed in tetrad basis (columns M = {t_E, r\:0302, \[Theta]\:0302, \[CurlyPhi]\:0302}) *)
Clear[SeedFromZ];
SeedFromZ[z_?NumericQ] := Module[{\[Psi] = psiFromZ[z], \[Alpha], \[Beta], P},
  \[Alpha] = 0.5 Sqrt[\[Psi]];  (* along r\:0302 *)
  \[Beta] = -Sqrt[\[Psi]];     (* along \[Theta]\:0302, \[CurlyPhi]\:0302 *)
  P = ConstantArray[0., {4, 4}];
  P[[2, 2]] = \[Alpha];    (* P^1_{r\:0302} *)
  P[[3, 3]] = \[Beta];    (* P^2_{\[Theta]\:0302} *)
  P[[4, 4]] = \[Beta];    (* P^3_{\[CurlyPhi]\:0302} *)
  P[[1, 1]] = 1/(\[Alpha] \[Beta]^2);
  Developer`ToPackedArray @ N @ P
];

$SQG1Loaded = True;



(* ::Package:: *)
(**)


(* ::Input:: *)
(* ========================== SQG2 . wl (Euclidean tetrad, metric-free) ========================== *)


(*tolerances/precision*)
If[!ValueQ[$SQGTol],$SQGTol=1.*^-12];
If[!ValueQ[$SQGWP],$SQGWP=80];

(*-----Pauli& internal basis T_A-----*)
pauli={{{0,1},{1,0}},{{0,-I},{I,0}},{{1,0},{0,-1}}};
tBasis[A_]:=If[A==1,I IdentityMatrix[2],pauli[[A-1]]];



(* ----- Euclidean mixed 't Hooft symbols (time slot = 1) ----- *)
(* A,B = 1..4 with 1 \[Congruent] t_E; spatial slots 2..4 map to a=1..3.
   \[Eta]^i_{ab} = \[CurlyEpsilon]_{iab},   \[Eta]^i_{a,1} = +\[Delta]_{ia},   \[Eta]^i_{1,a} = -\[Delta]_{ia}
   \.08ar\[Eta]^i_{ab} = \[CurlyEpsilon]_{iab}, \.08ar\[Eta]^i_{a,1} = -\[Delta]_{ia}, \.08ar\[Eta]^i_{1,a} = +\[Delta]_{ia} *)
Clear[EtaMixed, BarEtaMixed];
EtaMixed := EtaMixed = Module[{e3 = LeviCivitaTensor[3], \[Delta] = KroneckerDelta, E},
  E = ConstantArray[0., {3, 4, 4}];
  Do[E[[i, a+1, b+1]] = e3[[i, a, b]], {i, 1, 3}, {a, 1, 3}, {b, 1, 3}];
  Do[E[[i, a+1, 1]] =  \[Delta][i, a],       {i, 1, 3}, {a, 1, 3}];
  Do[E[[i, 1, a+1]] = -\[Delta][i, a],       {i, 1, 3}, {a, 1, 3}];
  Developer`ToPackedArray @ E
];
BarEtaMixed := BarEtaMixed = Module[{e3 = LeviCivitaTensor[3], \[Delta] = KroneckerDelta, E},
  E = ConstantArray[0., {3, 4, 4}];
  Do[E[[i, a+1, b+1]] = e3[[i, a, b]], {i, 1, 3}, {a, 1, 3}, {b, 1, 3}];
  Do[E[[i, a+1, 1]] = -\[Delta][i, a],       {i, 1, 3}, {a, 1, 3}];
  Do[E[[i, 1, a+1]] =  \[Delta][i, a],       {i, 1, 3}, {a, 1, 3}];
  Developer`ToPackedArray @ E
];


(*-----spacetime helpers-----*)
Pmu2x2[P_,\[Mu]_Integer?Positive]:=Sum[tBasis[A]*P[[A,\[Mu]]],{A,1,4}];

(*LOWER only (P stored with upper \[Mu]):P^A{}_\[Nu]=g_{\[Nu]\[Mu]} P^{A\[Mu]}*)

(*-----curvature F from P (for Urbantke)-----*)
Fhat[P_]:=Module[{F=ConstantArray[0.,{3,4,4}],Pm,Pn,comm},Do[Pm=Pmu2x2[P,\[Mu]];Pn=Pmu2x2[P,\[Nu]];comm=Pm . Pn-Pn . Pm;
Do[F[[i,\[Mu],\[Nu]]]=Im[Tr[pauli[[i]] . comm]];
F[[i,\[Nu],\[Mu]]]=-F[[i,\[Mu],\[Nu]]],{i,1,3}],{\[Mu],1,4},{\[Nu],\[Mu]+1,4}];
Developer`ToPackedArray[F]];

(* Clear[UrbantkeMetric];
UrbantkeMetric[P_, tol_:$SQGTol] := Module[{F, FF1, FF2, FF3, G, d, gDn},
  F = Fhat[P];                                (* F is 3\[Times]4\[Times]4 of antisymmetric 2-forms *)
  FF1 = F[[2]].F[[3]] - F[[3]].F[[2]];        (* \[CurlyEpsilon]_{23k} F^2 F^3, etc. *)
  FF2 = F[[3]].F[[1]] - F[[1]].F[[3]];
  FF3 = F[[1]].F[[2]] - F[[2]].F[[1]];
  G   = -(FF1.F[[1]] + FF2.F[[2]] + FF3.F[[3]])/6.;   (* \[Sqrt]g g = -(1/6) \[CurlyEpsilon] F F F *)
  G   = Chop[(G + Transpose[G])/2, tol];              (* symmetric numerically *)
  d   = Det @ N @ G;
  gDn = If[NumericQ[d] && d =!= 0, N @ Chop[G/Abs[d]^(1/4), tol], IdentityMatrix[4]];
  <|"densitized" -> G, "g" -> gDn, "ginv" -> Inverse[gDn]|>
]; *)

(*-----covariant internal Gram and A(\[CapitalOmega])-----*)
Clear[Sgram,AfromOmega,Skew];
Sgram[P_,tol_:$SQGTol]:= Developer`ToPackedArray@N@Chop[P . Transpose[P],tol];
AfromOmega[P_,Om_,tol_:$SQGTol]:= Developer`ToPackedArray@N@Chop[P . Om . Transpose[P],tol];
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
BuildSystem20[P_,S_,tol_:$SQGTol]:=Module[{Pnum=N@Chop[P,tol],\[CapitalPi]s,cols,Oe,He,dP,L18,trO,trH,L},
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
{Developer`ToPackedArray[L],\[CapitalPi]s}];

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
SolveOmegaH20[P_, S0_, f_, tol_ : $SQGTol, wp_ : $SQGWP] :=
 Module[{L,  \[CapitalPi]s, U, S, V, s, smax, d,  v, i0, zIdx, fUse},
  {L, \[CapitalPi]s} = BuildSystem20[P, S0, tol];

  {U, S, V} = SingularValueDecomposition[SetPrecision[L, wp]];
  s    = Diagonal[S];
  (* Print["after SVD: s=", N[s]]; *)
  (* pseudoinverse diag: nonzeros -> 1/s, zeros -> 0 *)
  smax = Max[s, 0.];
  d    = If[# > tol*smax, 1./#, 0.] & /@ s;
  (* indices of zero singulars: tail starting at first zero of d *)
  i0   = FirstPosition[d, 0., Missing["NotFound"]];
  zIdx = If[MissingQ[i0], {}, Range[i0[[1]], Length@d]];
  (* Print["null dim=", Length@zIdx]; *)
  (* match kick length to nullspace dim, or skip if none *)
  fUse = If[zIdx === {}, {}, Take[PadRight[Flatten@{f}, Length@zIdx, 0.], Length@zIdx]];
  If[fUse =!= {}, 
  v =  V[[All, zIdx]] . fUse;
  (* Print["v=", v]; *)
  ];

  <|
    "\[CapitalOmega]" -> ArrayReshape[v[[1 ;; 16]], {4, 4}],
    "h"               -> v[[17 ;; 20]],
    "\[CapitalPi]s"   -> \[CapitalPi]s
  |>
]


Clear[StepSDE20];
StepSDE20[P_,f_,d\[Theta]_:1.,stabEvery_:8,stepIndex_:1,tol_:$SQGTol,wp_:$SQGWP]:=
Module[{S,sol,Om,h,\[CapitalPi]s,Hint,Xi,dP,Pnew,det},
  S=Sgram[P,tol];
  (* Print["S=", MatrixForm@S]; *)
  sol=SolveOmegaH20[P,S,f,tol,wp];
  Om=sol["\[CapitalOmega]"];h=sol["h"];\[CapitalPi]s=sol["\[CapitalPi]s"];
  Hint=Sum[h[[k]] \[CapitalPi]s[[k]],{k,1,4}];
  (*internal H*)If[$AssertStepOnce&&stepIndex==1,AssertStepConsistency[P,Om,Hint,$SQGTol];$AssertStepOnce=False];
  Xi=XiFrom[P,S,Om,Hint,Automatic];
  (*internal 4\[Times]4*)dP=Om . P-P . Xi;
  (* Print["d\[Theta]dP=", MatrixForm@(d\[Theta] dP)]; *)
  (*mixed spacetime/internal*)
  Pnew=P+d\[Theta] dP;
  det=Det[Pnew];
  Assert[Abs[Im[det]] < 0.1 tol];
  Pnew=Pnew/Abs[det]^(1/4);
  (*enforce abelian row purely imaginary*)
  Pnew[[1,All]]= Re[Pnew[[1,All]]];
  (*optional:keep non-abelian rows real*)
  (*Pnew[[2;;4,All]]=Re[Pnew[[2;;4,All]]];*)
  <|"P"->Developer`ToPackedArray[Pnew],"\[CapitalOmega]"->Om,"h"->h|>
];

(*-----invariants (same API)-----*)
SafeSqrtSym[X_,tol_:$SQGTol]:=
Module[{S=Symmetrize[X],vals,vecs,vclip},
{vals,vecs}=Eigensystem[S];
(* Print["x_vals=",vals]; *)
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



(* ::Package:: *)
(* ====================== SQG3_euclid.wl ====================== *)
(* Gaussian Fourier forcing sampled along SDE trajectories.       *)
(* ============================================================ *)


$HistoryLength = 0; (* Avoid kernel history bloat. *)

If[!ValueQ[$SQGTol], $SQGTol = 1.*^-12];
If[!ValueQ[$SQGZ0],  $SQGZ0  = 1.0];

(* Gaussian Fourier profile generator (vector of length M). *)
Clear[MakeUProfile];
MakeUProfile[M_: 6, K_: 12, sigma_: 0.15, seed_: Automatic] :=
  BlockRandom[
    Module[{a, b, xi},
      a  = RandomVariate[NormalDistribution[0, sigma], {M, K}];
      b  = RandomVariate[NormalDistribution[0, sigma], {M, K}];
      xi = RandomVariate[NormalDistribution[0, sigma], M];
      Function[
        \[Theta],
        Table[
          xi[[m]] + Sum[a[[m, k]] Sin[k \[Theta]] + b[[m, k]] Cos[k \[Theta]], {k, 1, K}],
          {m, 1, M}
        ]
      ]
    ],
    RandomSeeding -> seed
  ];

(* Trajectory sampler: streams R-values out via user-supplied flush. *)
Clear[SampleWTrajectory];
SampleWTrajectory[
    z_?NumericQ,
    fFun_,
    s_: +1,
    nSteps_: 512,
    stabEvery_: 16,
    tol_: $SQGTol,
    thinStride_: 10,
    chunkSize_: 10^4
  ] := Module[{\[Theta]s, d\[Theta], P, f\[Theta], sol, k, R, stride = thinStride, chunk = chunkSize,
          buffer = {}, output = {}},

    If[!IntegerQ[stride] || stride <= 0, stride = 1];
    If[!IntegerQ[chunk]  || chunk  <= 0, chunk  = 10^4];

    \[Theta]s = Subdivide[0., 2 Pi, nSteps - 1];
    d\[Theta] = \[Theta]s[[2]] - \[Theta]s[[1]];
    P  = SeedFromZ[z];

    For[k = 1, k <= nSteps, k++,
      f\[Theta]  = fFun[\[Theta]s[[k]]];
      sol = StepSDE20[P, f\[Theta], d\[Theta], stabEvery, k, tol];
      P   = sol["P"];

      If[Mod[nSteps - k, stride] == 0,
        R = RfromP[P, z, s];
        AppendTo[buffer, R];
        If[Length[buffer] >= chunk,
          AppendTo[output, buffer];
          buffer = {};
        ];
      ];
    ];

    If[buffer =!= {}, AppendTo[output, buffer]];

    Developer`ToPackedArray @ Flatten[output]
  ];



Clear[LogKernelPool];
LogKernelPool[] := Module[{info},
  info = Quiet[
    ParallelEvaluate[{ $KernelID, $MachineName, $MachineAddresses, $ProcessID },
      DistributedContexts -> None],
    Parallel`NoSubKernels::nosub
  ];
  If[info === $Failed, info = {}];
  SQGPrint["Worker kernels detected: ", Length[info]];
  Scan[
    SQGPrint["  id=", #[[1]], " host=", #[[2]], " addresses=", #[[3]], " pid=", #[[4]]] &,
    info
  ];
];

Clear[EnsureRemoteKernel];
EnsureRemoteKernel[] := Module[{remoteHost = "100.111.63.24", remoteSpec, info, hasRemote},
  remoteHost = "mac-a102";
remoteSpec = RemoteKernel[
  "amigdal@" <> remoteHost,
  "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null amigdal@" <> remoteHost <>
    " '/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -mathlink -noinit'"
];

  info = Quiet[
    ParallelEvaluate[{ $KernelID, $MachineName, $MachineAddresses, $ProcessID },
      DistributedContexts -> None],
    Parallel`NoSubKernels::nosub
  ];
  If[info === $Failed || info === {},
    SQGPrint["No worker kernels active; launching default local kernels..."];
    Quiet[LaunchKernels[], Parallel`NoSubKernels::nosub];
    info = Quiet[
      ParallelEvaluate[{ $KernelID, $MachineName, $MachineAddresses, $ProcessID },
        DistributedContexts -> None],
      Parallel`NoSubKernels::nosub
    ];
    If[info === $Failed, info = {}];
  ];
  hasRemote = AnyTrue[
    info,
    (#[[2]] === remoteHost) || MemberQ[#[[3]], remoteHost] &
  ];

  SQGPrint["Current worker addresses: ", If[info === {}, "none", info[[All, 3]]]];
  If[!hasRemote,
    SQGPrint["Launching remote kernel on ", remoteHost, "..."];
    Check[Quiet[LaunchKernels[remoteSpec], Parallel`RemoteKernels`RemoteLaunch::fail], $Failed];
    Pause[1];
    info = Quiet[
      ParallelEvaluate[{ $KernelID, $MachineName, $MachineAddresses, $ProcessID },
        DistributedContexts -> None],
      Parallel`NoSubKernels::nosub
    ];
    If[info === $Failed, info = {}];
    hasRemote = AnyTrue[info, MemberQ[#[[3]], remoteHost] &];
    SQGPrint["Worker addresses after launch: ", If[info === {}, "none", info[[All, 3]]]];
  ];
  If[!hasRemote,
    SQGPrint["Warning: remote kernel at ", remoteHost, " is still unavailable."];
  ];
  (*LogKernelPool[];*)
];

Clear[NewRunDir];
NewRunDir[] := FileNameJoin[{Directory[], "SQG_runs"}];

Clear[RunNewSimulation];
RunNewSimulation[
    z_?NumericQ,
    nSteps_Integer: 512,
    M_Integer: 6,
    K_Integer: 12,
    sigma_: 0.15,
    s_: 1,
    stabEvery_: 15,
    thinStride_: 10,
    tol_: $SQGTol,
    nCycles_Integer: 10^4
  ] := Module[{useeds, rundir, params, paramsFile, rThinFile, chunkLists, data},

  rundir = NewRunDir[];
  If[!DirectoryQ[rundir], CreateDirectory[rundir]];
  SQGPrint["Run dir: ", rundir];

  params = <|
    "z" -> z,
    "nSteps" -> nSteps,
    "M" -> M,
    "K" -> K,
    "sigma" -> sigma,
    "s" -> s,
    "stabEvery" -> stabEvery,
    "thinStride" -> thinStride,
    "tol" -> tol,
    "nCycles" -> nCycles
  |>;

  paramsFile = FileNameJoin[{rundir, "params.txt"}];
  Export[paramsFile, params];
  Export[FileNameJoin[{rundir, "date_begin.txt"}], Now];

  rThinFile = FileNameJoin[{rundir, "rthin.bin"}];
  If[FileExistsQ[rThinFile], DeleteFile[rThinFile]];

  useeds = Table[j * ($KernelCount*10) + i, {j, nCycles}, {i, $KernelCount}];

  DistributeDefinitions[MakeUProfile, SampleWTrajectory, nCycles];
  SQGPrint["[master] Ensuring remote kernel availability..."];
  (*EnsureRemoteKernel[];*)
  (*LogKernelPool[];*)
  SQGPrint["Processing on ", $KernelCount, " kernels, seeds with thinStride=", thinStride,
    ", nCycles=", nCycles];

  Map[(
      chunkLists =     
        ParallelMap[
          Function[seed,
            Module[{f, samples},
              f = MakeUProfile[M, K, sigma, seed];
              samples = SampleWTrajectory[z, f, s, nSteps, stabEvery, tol, thinStride, nCycles];
              SQGPrint["[Kernel ", $KernelID, " @ ", $MachineName, "] completed seed ", seed,
                " with ", Length[samples], " samples."];
              samples
            ]
          ],
          #,
          ProgressReporting -> True
        ];

      (* Save chunk lists *)
      data = Developer`ToPackedArray @ N @ Flatten[chunkLists];
      SQGPrint["[master] writing ", Length[data], " samples to ", rThinFile];
      If[Length[data] > 0,
        Module[{fh = OpenWrite[rThinFile, BinaryFormat -> True]},
          BinaryWrite[fh, data, "Real64"];
          Close[fh];
        ]
      ];

      )&, 
    useeds
  ];

  Export[FileNameJoin[{rundir, "date_end.txt"}], Now];
  SQGPrint["Finished run ", rundir];
  LogKernelPool[];
];



(* Demo run parameters (adjust as needed). *)
On[Assert];
(*Plus@@ParallelTable[i,{i,24}]*)
RunNewSimulation[1., 1024, 6, 12, 1.1, 1, 1, 10, $SQGTol, 128];



