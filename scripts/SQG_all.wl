(* ::Package:: *)

(* ======================================================= *)
(* 0. Euclidean mixed 't Hooft symbols (time slot = 1)     *)
(* ======================================================= *)
Clear[EtaMixed,BarEtaMixed];

EtaMixed=Module[{e3=LeviCivitaTensor[3],\[Delta]=KroneckerDelta,E,i,a,b},E=ConstantArray[0.,{3,4,4}];
For[i=1,i<=3,i++,For[a=1,a<=3,a++,For[b=1,b<=3,b++,E[[i,a+1,b+1]]=e3[[i,a,b]];];];];
For[i=1,i<=3,i++,For[a=1,a<=3,a++,E[[i,a+1,1]]=\[Delta][i,a];
E[[i,1,a+1]]=-\[Delta][i,a];];];
Developer`ToPackedArray@E];

BarEtaMixed=Module[{e3=LeviCivitaTensor[3],\[Delta]=KroneckerDelta,E,i,a,b},E=ConstantArray[0.,{3,4,4}];
For[i=1,i<=3,i++,For[a=1,a<=3,a++,For[b=1,b<=3,b++,E[[i,a+1,b+1]]=e3[[i,a,b]];];];];
For[i=1,i<=3,i++,For[a=1,a<=3,a++,E[[i,a+1,1]]=-\[Delta][i,a];
E[[i,1,a+1]]=\[Delta][i,a];];];
Developer`ToPackedArray@E];

Aeta := Table[EtaMixed[[a, ;;, ;;]], {a, 1, 3}];

(* ======================================================= *)
(* 1. Seed P(z)                                            *)
(* ======================================================= *)
Clear[psiFromZ, SeedFromZ];
psiFromZ[z_?NumericQ] := 1./(2. z);

SeedFromZ[z_?NumericQ] := Module[{\[Psi] = psiFromZ[z], \[Alpha], \[Beta], P},
  \[Alpha] = 0.5 Sqrt[\[Psi]]; \[Beta] = -Sqrt[\[Psi]];
  P = ConstantArray[0., {4, 4}];
  P[[2, 2]] = \[Alpha]; P[[3, 3]] = \[Beta]; P[[4, 4]] = \[Beta];
  P[[1, 1]] = 1/(\[Alpha] \[Beta]^2);
  Developer`ToPackedArray @ N @ P
];

Clear[Skew];
Skew[X_] := (X - Transpose[X])/2;

(* ======================================================= *)
(* 2. f(P) and \[Delta]f(P,\[Delta]P)                                   *)
(* ======================================================= *)
Clear[fMatrix, dfMatrixByDeltaP];

fMatrix[P_] := Module[{Amats = Aeta},
  Table[
    1/2 Sum[
      LeviCivitaTensor[3][[i, j, k]] *
        ( P[[j, ;;]] . Amats[[a]] . P[[;;, k]] ),
      {j, 1, 3}, {k, 1, 3}
    ],
    {i, 1, 3}, {a, 1, 3}
  ]
];

dfMatrixByDeltaP[P_, dP_] := Module[{Amats = Aeta},
  Table[
    1/2 Sum[
      LeviCivitaTensor[3][[i, j, k]] *
        ( dP[[j, ;;]] . Amats[[a]] . P[[;;, k]] + 
          P[[j, ;;]] . Amats[[a]] . dP[[;;, k]] ),
      {j, 1, 3}, {k, 1, 3}
    ],
    {i, 1, 3}, {a, 1, 3}
  ]
];

(* ======================================================= *)
(* 3. Stable Xi(P,\[CapitalOmega],H) via clipped eigen-inverse           *)
(* ======================================================= *)
Clear[XiFromOmegaHClip];
XiFromOmegaHClip[P_, \[CapitalOmega]_, H_, clip_:1.*^-12] := Module[
  {S, rhs, vals, vecs, \[Lambda]max, inv\[Lambda], invS},
  S   = Transpose[P] . P;
  rhs = Skew[Transpose[P] . \[CapitalOmega] . P] + H;
  {vals, vecs} = Eigensystem[S];
  \[Lambda]max = Max[Abs[vals]];
  inv\[Lambda] = If[# < clip*\[Lambda]max, 0., 1./#] & /@ vals;
  invS = vecs . DiagonalMatrix[inv\[Lambda]] . Transpose[vecs];
  - invS . rhs
];

(* ======================================================= *)
(* 4. Duality constraints (18) + extra 4 rows              *)
(* ======================================================= *)
Clear[DualityRowsForDeltaP, ExtraRowsFromOmega];

DualityRowsForDeltaP[P_, dP_] := Module[{Amats = Aeta, rows1, rows2},
  rows1 = Flatten @ Table[
    Sum[LeviCivitaTensor[3][[i, j, k]] * ( P[[j, ;;]] . Amats[[a]] . dP[[;;, k]] ),
      {j, 1, 3}, {k, 1, 3}],
    {a, 1, 3}, {i, 1, 3}
  ];
  rows2 = Flatten @ Table[
    Sum[LeviCivitaTensor[3][[i, j, k]] * ( dP[[j, ;;]] . Amats[[a]] . P[[;;, k]] ),
      {j, 1, 3}, {k, 1, 3}],
    {a, 1, 3}, {i, 1, 3}
  ];
  Join[rows1, rows2]
];

ExtraRowsFromOmega[\[CapitalOmega]_] := Join[{Tr[\[CapitalOmega]]},
  {\[CapitalOmega][[2,3]]-\[CapitalOmega][[3,2]], \[CapitalOmega][[2,4]]-\[CapitalOmega][[4,2]], \[CapitalOmega][[3,4]]-\[CapitalOmega][[4,3]]}];

(* ======================================================= *)
(* 5. Parameterization (26 = 16 \[CapitalOmega] + 10 H)                  *)
(* ======================================================= *)
Clear[upperPairs, unvecOmega, unpackH, ParamToOmegaH];
upperPairs = {{1,1},{1,2},{1,3},{1,4},{2,2},{2,3},{2,4},{3,3},{3,4},{4,4}};

unvecOmega[\[Omega]_List] /; Length[\[Omega]]==16 := Partition[\[Omega], 4];

unpackH[h_List] /; Length[h]==10 := Module[{H = ConstantArray[0., {4, 4}]},
  Do[
    With[{\[Mu] = upperPairs[[k, 1]], \[Nu] = upperPairs[[k, 2]]},
      H[[\[Mu], \[Nu]]] = h[[k]]; H[[\[Nu], \[Mu]]] = h[[k]];
    ], {k, 10}];
  H
];

ParamToOmegaH[g_List] /; Length[g]==26 := Module[{\[Omega] = Take[g, 16], h = Take[g, -10]},
  {unvecOmega[\[Omega]], unpackH[h]}
];

(* ======================================================= *)
(* 6. Build E(P)                                           *)
(* ======================================================= *)
Clear[DeltaPFromParam, BuildEParam];

DeltaPFromParam[P_, j_Integer, clip_:1.*^-12] := Module[{g = ConstantArray[0., 26], \[CapitalOmega], H},
  g[[j]] = 1.;
  {\[CapitalOmega], H} = ParamToOmegaH[g];
  \[CapitalOmega] . P - P . XiFromOmegaHClip[P, \[CapitalOmega], H, clip]
];

BuildEParam[P_, clip_:1.*^-12] := Module[{cols},
  cols = Table[
    Module[{dP = DeltaPFromParam[P, j, clip], \[CapitalOmega], H, base, extra},
      {\[CapitalOmega], H} = ParamToOmegaH[UnitVector[26, j]];
      base  = DualityRowsForDeltaP[P, dP];
      extra = ExtraRowsFromOmega[\[CapitalOmega]];
      Join[base, extra]
    ], {j, 26}];
  Transpose @ Developer`ToPackedArray @ cols
];

(* ======================================================= *)
(* 7. Nullspace report                                     *)
(* ======================================================= *)
Clear[NullspaceReport];
NullspaceReport[P_, tol_:1.*^-10, clip_:1.*^-12] := Module[
  {E = BuildEParam[P, clip], U, S, V, s, eps, zidx},
  {U, S, V} = SingularValueDecomposition[N @ E];
  s   = Diagonal[S]; eps = tol * Max[Join[s, {0.}]];
  zidx = Flatten @ Position[s, _?(# <= eps &)];
  <|"E"->E, "singulars"->s, "tol"->eps, "Z"->Length[zidx], "Vnull"->V[[All, zidx]]|>
];

(* ======================================================= *)
(* 8. R(P) and dR(P;\[Delta]P)                                    *)
(* ======================================================= *)
Clear[RofP, dRLinear];
RofP[P_] := 4 * Total @ SingularValueList[N @ fMatrix[P]];

dRLinear[P_, dP_] := Module[{f = fMatrix[P], df = dfMatrixByDeltaP[P, dP], U, S, V, Rpol},
  {U, S, V} = SingularValueDecomposition[N @ f];
  Rpol = U . Transpose[V];  (* (ff^T)^(-1/2) f *)
  8 * Tr[Rpol . Transpose[df]]
];

(* ======================================================= *)
(* 9. OneStep                                              *)
(* ======================================================= *)
Clear[OneStep];
OneStep[P_, G_, tol_:1.*^-10, clip_:1.*^-12] := Module[
  {rep, Vz, Z, len, Gadj, g, \[CapitalOmega], H, Xi, dP, Pnew, Rnew, dR},
  rep = NullspaceReport[P, tol, clip];
  If[rep["Z"] == 0, Return[{P, Missing["NoNull"], Missing["NoNull"]}]];
  Vz  = rep["Vnull"];  Z = Dimensions[Vz][[2]];  len = Length[G];
  Gadj = If[len === Z, G,
            If[len > Z, Take[G, Z], Join[G, ConstantArray[0., Z - len]]]];
  g     = Vz . Gadj;
  {\[CapitalOmega], H} = ParamToOmegaH[g];
  Xi    = XiFromOmegaHClip[P, \[CapitalOmega], H, clip];
  dP    = \[CapitalOmega] . P - P . Xi;
  Pnew  = P + dP;
  Rnew  = RofP[Pnew];
  dR    = dRLinear[P, dP];
  {Pnew, Rnew, dR}
];

(* ======================================================= *)
(* RDerivative[P,G]: analytic derivative along direction G *)
(* ======================================================= *)

Clear[RDerivative];
RDerivative[P_, G_, tol_:1.*^-10, clip_:1.*^-12] := Module[
  {rep, Vz, Z, len, Gadj, g, \[CapitalOmega], H, Xi, dP, dR},
  rep = NullspaceReport[P, tol, clip];
  If[rep["Z"] == 0, Return[Missing["NoNull"]]];
  Vz  = rep["Vnull"];  Z = Dimensions[Vz][[2]];  len = Length[G];
  Gadj = If[len === Z, G,
            If[len > Z, Take[G, Z], Join[G, ConstantArray[0., Z - len]]]];
  g     = Vz . Gadj;
  {\[CapitalOmega], H} = ParamToOmegaH[g];
  Xi    = XiFromOmegaHClip[P, \[CapitalOmega], H, clip];
  dP    = \[CapitalOmega] . P - P . Xi;
  dR    = dRLinear[P, dP];
  dR
];


Dimensions[EtaMixed]



z = 1.;
P0 = SeedFromZ[z];


(* 1. Is P0 numeric? *)
MatrixQ[P0, NumericQ]


(* 2. Check XiFromOmegaHClip directly on P0 with a simple \[CapitalOmega],H *)
\[CapitalOmega]test = ConstantArray[0., {4,4}]; Htest = ConstantArray[0., {4,4}];
XiFromOmegaHClip[P0, \[CapitalOmega]test, Htest] // MatrixForm


(* 3. BuildEParam should give a 22\[Times]26 numeric matrix *)
E = BuildEParam[P0];
Dimensions[E]
MatrixQ[E, NumericQ]


(* 4. SingularValueDecomposition directly *)
SingularValueDecomposition[E][[2]] // Diagonal



z = 1.;
P0 = SeedFromZ[z];
rep = NullspaceReport[P0, 1.*^-10];
{rep["Z"], rep["singulars"]}


Z   = rep["Z"];
G   = 10^-6 * Table[(-1)^k, {k, Z}];

{P1, R1, dR} = OneStep[P0, G];
{Det[P0] // N, Det[P1] // N, RofP[P0], R1, dR}



(* ::Package:: *)
(* ====================== SQG3_euclid.wl ====================== *)
(* Gaussian Fourier forcing sampled along SDE trajectories.       *)
(* ============================================================ *)


$HistoryLength = 0; (* Avoid kernel history bloat. *)

If[!ValueQ[$SQGTol], $SQGTol = 1.*^-12];
If[!ValueQ[$SQGZ0],  $SQGZ0  = 1.0];

(* Gaussian Fourier profile generator (vector of length M). *)
Clear[MakeUProfile];
MakeUProfile[M_: 3, K_: 12, sigma_: 0.15, seed_: Automatic] :=
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
      f\[Theta]  = d\[Theta]*fFun[\[Theta]s[[k]]];
      {P,R} = OneStep[P,f\[Theta]];
      Print["k=",k,"R=,",R];
      (*P = RK4StepCovariant[P, \[Theta]s[[k]], fFun, d\[Theta]];*)
      If[Mod[nSteps - k, stride] == 0,
        AppendTo[buffer, R];
        If[Length[buffer] >= chunk,
          AppendTo[output, buffer];
          SQGPrint["mean R , err =",Mean[buffer], StandardDeviation[buffer]]; 
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
NewRunDir[n_] := FileNameJoin[{Directory[], "SQG_runs_"<>ToString[n]}];

Clear[RunNewSimulation];
RunNewSimulation[
    n_?IntegerQ,
    nSteps_Integer: 512,
    M_Integer: 3,
    K_Integer: 12,
    sigma_: 0.15,
    s_: 1,
    stabEvery_: 15,
    thinStride_: 10,
    tol_: $SQGTol,
    nCycles_Integer: 10^4
  ] := Module[{useeds, rundir, params, paramsFile, rThinFile, chunkLists, data, fhThin, totalSamples = 0},

  rundir = NewRunDir[n];
  z = 2.^(-n);
  If[!DirectoryQ[rundir], CreateDirectory[rundir]];
  Print["Run dir: ", rundir];

  params = <|
    "n" -> n,
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
  numSeeds = $KernelCount;
  useeds = Table[j * ($KernelCount*10) + i, {j, nCycles}, {i, numSeeds}];

  DistributeDefinitions[MakeUProfile, SampleWTrajectory, nCycles];
  SQGPrint["[master] Ensuring remote kernel availability..."];
  (*EnsureRemoteKernel[];*)
  (*LogKernelPool[];*)
  SQGPrint["Processing on ", $KernelCount, " kernels, seeds with thinStride=", thinStride,
    ", nCycles=", nCycles];

  fhThin = OpenAppend[rThinFile, BinaryFormat -> True];
  If[fhThin === $Failed,
    SQGPrint["Error: unable to open ", rThinFile, " for appending."];
    Return[$Failed];
  ];

  Map[(
      chunkLists =     
        ParallelMap[
          Function[seed,
            Module[{f, samples},
              f = MakeUProfile[M, K, sigma, seed];
              SampleWTrajectory[z, f, s, nSteps, stabEvery, tol, thinStride, nCycles]
            ]
          ],
          #,
          ProgressReporting -> True
        ];

      (* Save chunk lists *)
      data = Developer`ToPackedArray @ N @ Flatten[chunkLists];
      If[Length[data] > 0,
        BinaryWrite[fhThin, data, "Real64"];
        Flush[fhThin];
        totalSamples += Length[data];
        Print["[master] appended ", Length[data], " samples (total ", totalSamples,"; ", Mean[data],"+/-", StandardDeviation[data]," to ", rThinFile];
      ];

      )&, 
    useeds
  ];

  Close[fhThin];

  Export[FileNameJoin[{rundir, "date_end.txt"}], Now];
  SQGPrint["Finished run ", rundir];
  LogKernelPool[];
];



(* 0) Clean slate on master *)
$HistoryLength = 0;
ClearAll[$LRCache, BuildSymbolicLR, StageParamJacobian,
  NullspaceParam, UnpackParamStep, DeltaXiFrom, RK4StepCovariant,
  Sproj4, Sgram, Skew, P2x2FromColumn, dP2x2FromColumn];

(* 1) (Re)define pauli, BarEtaMixed, helpers, and all functions *)
(* ... your definitions here ... *)

(* 2) Build the symbolic Grad cache ON MASTER only *)
BuildSymbolicLR[];

(* 3) Ensure kernels are up *)
If[$KernelCount == 0, LaunchKernels[]];

(* 4) Distribute ONLY function definitions & small symbols *)
DistributeDefinitions[
  pauli, BarEtaMixed,
  Skew, Sgram, Sproj4,
  P2x2FromColumn, dP2x2FromColumn,
  BuildSymbolicLR, StageParamJacobian, NullspaceParam,
  UnpackParamStep, DeltaXiFrom, RK4StepCovariant
];

(* 5) Build the cache on each worker (DON\[CloseCurlyQuote]T ship $LRCache) *)
ParallelEvaluate[BuildSymbolicLR[]];

(* 6) Now run your parallel simulation *)
Do[RunNewSimulation[n, 2048, 3, 64, 0.25, 1, 1, 10, $SQGTol, 128], {n, 0, 1}];




