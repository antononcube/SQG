(* ::Package:: *)

 Clear[Skew]; 
  Skew[X_] := (X - Transpose[X])/2;
  
  Clear[Ep]; 
  Ep = Table[Table[Which[j <= 3 && k <= 3, LeviCivitaTensor[3][[i,j,k]], 
       j == 4 && k <= 3, KroneckerDelta[i, k], j <= 3 && k == 4, 
       KroneckerDelta[i, j], True, 0], {j, 4}, {k, 4}], {i, 3}]; 
  Clear[EtaMixed, BarEtaMixed]; 
  EtaMixed = Table[Table[Which[j <= 3 && k <= 3, LeviCivitaTensor[3][[i,j,
        k]], j <= 3 && k == 4, KroneckerDelta[i, j], j == 4 && k <= 3, 
       -KroneckerDelta[i, k], True, 0], {j, 4}, {k, 4}], {i, 3}]; 
  BarEtaMixed = Table[Table[Which[j <= 3 && k <= 3, LeviCivitaTensor[3][[i,j,
        k]], j <= 3 && k == 4, -KroneckerDelta[i, j], j == 4 && k <= 3, 
       KroneckerDelta[i, k], True, 0], {j, 4}, {k, 4}], {i, 3}]; 
  Aeta := Table[EtaMixed[[a]], {a, 1, 3}]; 
  BAeta := Table[BarEtaMixed[[a]], {a, 1, 3}]; 
  Clear[BuildT]; 
  BuildT[] := Module[{P, Q, dualityeqs, traceeqs, sk, eqs, varsP, varsQ, 
      coeffs, Mlin, T}, P = Array[p, {4, 4}]; Q = Array[q, {4, 4}]; 
      dualityeqs = Flatten[Table[Tr[Ep[[i]] . P . BAeta[[a]] . Transpose[Q]], 
         {i, 3}, {a, 3}]]; sk = Skew[P . Transpose[Q]]; 
      traceeqs = {sk[[1,2]], sk[[1,3]], sk[[1,4]], sk[[2,3]], sk[[2,4]], 
        sk[[3,4]]}; eqs = Join[dualityeqs, traceeqs]; varsP = Flatten[P]; 
      varsQ = Flatten[Q]; coeffs = CoefficientArrays[eqs, varsQ]; 
      Mlin = coeffs[[2]]; T = Table[Coefficient[Mlin[[row,col]], varsP[[d]]], 
        {row, 15}, {col, 16}, {d, 16}]; T]; 
        Tcoeff = BuildT[]; 
  Clear[BuildM]; 
  BuildM[P_] := Module[{pflat = Flatten[P]}, 
     Table[pflat . Tcoeff[[row,col]], {row, 15}, {col, 16}]];

(* Canonical nullspace: always 16 x k (columns are null vectors) *)
Clear[NullspaceV];
NullspaceV[P_, tol_: 1.*^-10] := Module[{M, ns},
  M  = N @ BuildM[P];                      (* 15 x 16 *)
  ns = NullSpace[M, Tolerance -> tol];     (* either {16} or {k,16} *)

  Which[
    ns === {}, Return[$Failed],
    VectorQ[ns, NumericQ],     (* {16} case: a single 16-vector *)
      Transpose @ {ns},        (* \[RightArrow] 16 x 1 *)
    MatrixQ[ns, NumericQ],     (* {k,16} case: list-of-vectors *)
      Transpose[ns],           (* \[RightArrow] 16 x k *)
    True,
      Return[$Failed]
  ]
];

Clear[QfromP];
QfromP[P_, tol_: 1.*^-10] := Module[{Vnull, k, G, qflat, Q, nP, nQ},
  Vnull = NullspaceV[P, tol];
  If[Vnull === $Failed, Return[Null]];

  (* Now Vnull is 16 x k by construction *)
  k = Dimensions[Vnull][[2]];
  (* Print["[QfromP] Nullity = ", k]; *)

  (* random combo of the k modes (k==1 \[RightArrow] scalar) *)
  G = RandomVariate[NormalDistribution[0,1], k];

  qflat = Vnull . G;                    (* length 16 *)
  If[Length[qflat] =!= 16, Return[Null]];

  Q  = ArrayReshape[qflat, {4,4}];
  nQ = Tr[Q . Transpose[Q]]; If[nQ == 0., Return[Null]];
  nP = Tr[P . Transpose[P]];

  Sqrt[nP/nQ] * Q
];




Clear[ fMatrix ];
fMatrix[P_, Q_] := 1/4Table[Tr[Ep[[i]] . P . Aeta[[a]] . Transpose[Q]], 
      {i, 3}, {a, 3}];


(* ======================= InitPQ Core Routines ======================= *)

Clear[SPQfTarget];
SPQfTarget[z_?NumericQ] := Module[{psi = 1./(2. z)},
  {{2 psi, 0., 0.}, {0., -psi, 0.}, {0., 0., -psi}}
];

Clear[SPQBuildLR];
SPQBuildLR[] := Module[{Llist = {}, Rlist = {}, I4 = IdentityMatrix[4], pairs, Eij},
  Do[
    AppendTo[Llist, Normal @ Ep[[i]]];
    AppendTo[Rlist, Normal @ BAeta[[a]]],
    {i, 1, 3}, {a, 1, 3}
  ];
  pairs = {{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}};
  Do[
    Eij = SparseArray[{{pairs[[k, 1]], pairs[[k, 2]]} -> 1.}, {4, 4}];
    AppendTo[Llist, Normal @ 0.5 (Eij - Transpose[Eij])];
    AppendTo[Rlist, I4],
    {k, Length[pairs]}
  ];
  Do[
    AppendTo[Llist, Normal @ Ep[[i]]];
    AppendTo[Rlist, Normal @ Aeta[[a]]],
    {i, 1, 3}, {a, 1, 3}
  ];
  <|"L" -> Llist, "R" -> Rlist|>
];

Clear[SPQBuildB];
SPQBuildB[z_?NumericQ] := Module[{F = SPQfTarget[z]},
  N @ Join[ConstantArray[0., 9], ConstantArray[0., 6], Flatten[4*F]]
];

Clear[SPQAQ];
SPQAQ[P_, L_, R_, b_] := Module[{rows, ok},
  rows = MapThread[
    Developer`ToPackedArray @ Flatten @ Normal[#1 . P . #2] &,
    {L, R}
  ];
  ok = VectorQ[rows, VectorQ[#, NumericQ] && Length[#] == 16 &];
  If[!ok,
    Print["[SPQAQ] malformed rows: row count=", Length[rows],
      ", unique lengths=", DeleteDuplicates[Length /@ rows]];
    Print["  example row lengths=", Take[Length /@ rows, UpTo[5]]];
    Print["  first problematic row=", First[Select[rows, Length[#] =!= 16 &]]];
    Abort[];
  ];
  {rows, Developer`ToPackedArray @ N @ b}
];

Clear[SPQAP];
SPQAP[Q_, L_, R_, b_] := Module[{rows, ok},
  rows = MapThread[
    Developer`ToPackedArray @ Flatten @ Normal[Transpose[#2 . Transpose[Q] . #1]] &,
    {L, R}
  ];
  ok = VectorQ[rows, VectorQ[#, NumericQ] && Length[#] == 16 &];
  If[!ok,
    Print["[SPQAP] malformed rows: row count=", Length[rows],
      ", unique lengths=", DeleteDuplicates[Length /@ rows]];
    Print["  example row lengths=", Take[Length /@ rows, UpTo[5]]];
    Print["  first problematic row=", First[Select[rows, Length[#] =!= 16 &]]];
    With[{sampleL = First[L], sampleR = First[R]},
      Print["  dims sample: L -> ", Dimensions[sampleL],
        " R -> ", Dimensions[sampleR]]
    ];
    Abort[];
  ];
  {rows, Developer`ToPackedArray @ N @ b}
];

Clear[SPQPinRowsP];
SPQPinRowsP[p4_List] /; Length[p4] == 4 := Module[{rows = {}, rhs = {}, E},
  Do[
    E = SparseArray[{{4, j} -> 1.}, {4, 4}];
    AppendTo[rows, Developer`ToPackedArray @ Flatten @ Normal[E]];
    AppendTo[rhs, p4[[j]]],
    {j, 1, 4}
  ];
  {rows, rhs}
];

Clear[SPQSolveQ, SPQSolveP];
SPQSolveQ[P_, LR_Association, b_] := Module[{L = LR["L"], R = LR["R"], A, B, q},
  {A, B} = SPQAQ[P, L, R, b];
  If[!MatrixQ[A],
    Print["[SPQSolveQ] non-rectangular matrix A; dims=", Dimensions[A]];
    Abort[];
  ];
  q = LeastSquares[A // Normal, B, Method -> "Direct"];
  ArrayReshape[q, {4, 4}]
];

SPQSolveP[Q_, p4_List, LR_Association, b_] := Module[{L = LR["L"], R = LR["R"], A, B, pinR, pinB, p},
  {A, B} = SPQAP[Q, L, R, b];
  {pinR, pinB} = SPQPinRowsP[p4];
  If[!MatrixQ[A],
    Print["[SPQSolveP] non-rectangular matrix A; dims=", Dimensions[A],
      ", pin dims=", Dimensions[pinR]];
    Abort[];
  ];
  p = LeastSquares[Normal@Join[A, pinR], Join[B, pinB], Method -> "Direct"];
  ArrayReshape[p, {4, 4}]
];

Clear[SPQSolvePQ];
SPQSolvePQ[z_?NumericQ, LR_Association, maxIter_: 60, tol_: 1.*^-10] :=
  Module[{b = SPQBuildB[z], P, p4, Q, res, it = 0},
    p4 = Normalize@RandomVariate[NormalDistribution[0, 1], 4];
    P = ConstantArray[0., {4, 4}];
    P[[4]] = p4;
    Q = ConstantArray[0., {4, 4}];
    While[it < maxIter,
      it++;
      Q = SPQSolveQ[P, LR, b];
      P = SPQSolveP[Q, p4, LR, b];
      res = Norm @ (MapThread[Tr[#1 . P . #2 . Transpose[Q]] &, {LR["L"], LR["R"]}] - b);
      If[res <= tol, Break[]];
    ];
    <|"P" -> P, "Q" -> Q, "P4" -> p4, "Iterations" -> it, "Residual" -> res|>
  ];


InitParallel[] := (
  LaunchKernels[];
  $TcoeffMaster = BuildT[];
  DistributeDefinitions[Ep, EtaMixed, BarEtaMixed, Aeta, BAeta,
    BuildT, BuildM, QfromP, fMatrix, SeedFromZ, Skew, $TcoeffMaster];
  ParallelEvaluate[Tcoeff = $TcoeffMaster];
);
InitParallel[];

(* ======================= Example Usage ======================= *)
LR = SPQBuildLR[];
DistributeDefinitions[LR];
  Clear[OneRun]; 
  OneRun[z_, numdata_, tol_:1.*^-16] := 
    Module[{P, sol, RRR = {}, Q, R}, sol = SPQSolvePQ[z, LR, 100, tol]; 
      P = sol["Q"]; Do[Q = QfromP[P, tol]; If[Q != {}, 
         R = 4*Total[SingularValueList[N[fMatrix[P, Q]]]]; AppendTo[RRR, R]; 
          P = Q, sol = SPQSolvePQ[z, LR, 100, tol]; P = sol["Q"]], 
       {numdata}]; RRR];
RunBatch[z_:1., cycles_:16, samples_:1000, tol_:1.*^-16] :=
  Module[{results},
    results = ParallelTable[OneRun[z, samples, tol], {cycles}];
    Sort[Flatten[results]]
  ];


ClearAll[PlotTailNormalized, MakeTailsNormalized, normBy];

(*normalization helper*)
normBy[v_List, "Median"] := Median[v];
normBy[v_List, "Mean"] := Mean[v];

(*----------single dataset:{R,z}----------*)
Options[PlotTailNormalized] = {"NormalizeBy" -> 
    "Mean"};   (*"Mean"|"Median"*)

PlotTailNormalized[{R_, z_}, OptionsPattern[]] := 
  Module[{x, n, \[Mu], tail, lbl, how = OptionValue["NormalizeBy"]}, 
   x = Select[Flatten@Normal@N@{R}, NumberQ];
   If[Length[x] < 2, 
    Return@Graphics[{}, Frame -> True, 
      FrameLabel -> {"R/\[Mu]", "1 - CDF (log)"}]];
   n = Length[x];
   \[Mu] = normBy[x, how];
   tail = (n - Range[n] + 1.)/(n + 1.);(*Weibull:never 0*)
   lbl = "z=" <> ToString@NumberForm[z, {8, 3}];
   ListPlot[Transpose[{x/\[Mu], tail}], Joined -> True, 
    Frame -> True, 
    FrameLabel -> {"R / " <> ToString[how], "1 - CDF (log)"}, 
    ScalingFunctions -> {"Linear", "Log"}, PlotRange -> All, 
     PlotLabel->"Weyl curvature 1-CDF",
    PlotLegends -> {lbl}, ImageSize -> Large]];

(*----------multiple datasets:{R1,z1},{R2,z2},...----------*)
Options[MakeTailsNormalized] = Options[PlotTailNormalized];

MakeTailsNormalized[pairs__List, OptionsPattern[]] := 
  Module[{items = {pairs}, curves = {}, labels = {}, 
    how = OptionValue["NormalizeBy"]}, 
   Do[Module[{R = items[[i, 1]], z = items[[i, 2]], x, n, \[Mu]}, 
     x = Select[Flatten@Normal@N@{R}, NumberQ];
     If[Length[x] >= 2, n = Length[x]; \[Mu] = normBy[x, how];
      AppendTo[curves, 
       Transpose[{x/\[Mu], (n - Range[n] + 1.)/(n + 1.)}]];
      AppendTo[labels, 
       "z=" <> ToString@NumberForm[z, {8, 3}]<>", <R>=" <> ToString@NumberForm[\[Mu], {8, 3}]]
       ]], {i, Length[items]}];
   If[curves === {}, 
    Return@Graphics[{}, Frame -> True, 
      FrameLabel -> {"R / " <> ToString[how], "1 - CDF (log)"}]];
   ListPlot[curves, Joined -> True, Frame -> True, 
    FrameLabel -> {"R / " <> ToString[how], "1 - CDF (log)"}, 
    ScalingFunctions -> {"Linear", "Log"}, PlotRange -> All, 
    PlotLabel->"Weyl curvature 1-CDF",
    PlotLegends -> labels, ImageSize -> Large]];
    
r=.;



(* helper: finite numeric predicate *)
isFiniteQ[expr_] := NumericQ[expr] && FreeQ[expr, _DirectedInfinity | Indeterminate];

Clear[FitMuAndGammas];

isFiniteQ[expr_] := NumericQ[expr] && FreeQ[expr, _DirectedInfinity | Indeterminate];

(* Input: curves = list (or list-of-lists) of quads {z, m, M, mu} *)
FitMuAndGammas[curves_] := Module[
  {flat, quads, rows, mudata, mudataFinite, x,
   mumodel, muindex, muSE,
   grouped, uniqueM, arrays, fits, validM, gamma, gammaSE,
   tau, tauSE, plotGamma, plotTau, table},

  (* Flatten and keep numeric quads *)
  flat  = Flatten[curves, 1];
  quads = Cases[flat, {z_?NumericQ, m_?NumericQ, M_?NumericQ, mu_?NumericQ}];
  If[quads === {}, Return[$Failed, Module]];

  (* Mean index: fit Log[mu] vs -Log[z] *)
  mudata = ({ -Log[#[[1]]], Log[#[[4]]]} & /@ quads);
  mudataFinite = Select[mudata, isFiniteQ[#[[1]]] && isFiniteQ[#[[2]]] &];
  If[mudataFinite === {},
    muindex = Missing["InsufficientData"]; muSE = Missing["NA"],
    mumodel = Quiet @ LinearModelFit[mudataFinite, {1, x}, x];
    muindex = mumodel["BestFitParameters"][[2]];
    muSE    = mumodel["ParameterErrors"][[2]];
  ];

  (* Build {m, x=-Log[z], y=Log[M]} and keep finite *)
  rows = Reap[
    Do[
      With[{m = quads[[i,2]], xx = -Log[quads[[i,1]]], yy = Log[quads[[i,3]]]},
        If[ isFiniteQ[xx] && isFiniteQ[yy], Sow[{m, xx, yy}] ]
      ],
      {i, Length[quads]}
    ]
  ][[2,1]];
  If[rows === {}, Return[$Failed, Module]];

  (* Group by m \[RightArrow] {{x,y},...} *)
  grouped = KeySort @ GroupBy[rows, First -> ({#[[2]], #[[3]]} &)];
  uniqueM = Keys[grouped]; arrays = Values[grouped];

  (* Fit y = a + gamma x for each m *)
  fits = AssociationThread[
    uniqueM ->
      Map[
        If[Length[#] >= 2 && Variance[#[[All,1]]] > 0,
           Quiet @ LinearModelFit[#, {1, x}, x],
           Missing["InsufficientData"]
        ] &,
        arrays
      ]
  ];

  validM = Select[uniqueM, MatchQ[fits[#], _LinearModelFit] &];

  gamma   = AssociationThread[validM, (fits[#]["BestFitParameters"][[2]] & /@ validM)];
  gammaSE = AssociationThread[validM, (fits[#]["ParameterErrors"][[2]]    & /@ validM)];

  (* Combined exponent tau(m) = m*muindex + gamma(m), with propagated SE *)
  tau  = AssociationThread[validM, (m |-> (m*muindex + gamma[m])) /@ validM];
  tauSE = AssociationThread[validM, (m |-> Sqrt[(m*muSE)^2 + (gammaSE[m])^2]) /@ validM];

  (* Plots *)
  plotGamma = ListPlot[
    ({#, Around[gamma[#],  gammaSE[#]]} & /@ Sort[validM]),
    Frame -> True, FrameLabel -> {"m", "\!\(\*SubscriptBox[\(\[Gamma]\), \(m\)]\)"},
    PlotMarkers -> Automatic, PlotTheme -> "Scientific", ImageSize -> 520
  ];

  plotTau = ListPlot[
    ({#, Around[tau[#], tauSE[#]]} & /@ Sort[validM]),
    Frame -> True, FrameLabel -> {"m", "\!\(\*SubscriptBox[\(\[Tau]\), \(m\)]\) = m \[Mu] + \[Gamma]"},
    PlotMarkers -> Automatic, PlotTheme -> "Scientific", ImageSize -> 520
  ];

  table = Grid[
    Prepend[
      ({#, gamma[#], gammaSE[#], tau[#], tauSE[#], fits[#]["RSquared"]} & /@ Sort[validM]),
      {"m", "\[Gamma](m)", "StdErr(\[Gamma])", "\[Tau](m)", "StdErr(\[Tau])", "R^2(\[Gamma]-fit)"}
    ],
    Frame -> All
  ];

  <|
    "muindex"  -> Around[muindex, muSE],
    "Gamma"    -> KeySort[gamma],
    "GammaErr" -> KeySort[gammaSE],
    "Tau"      -> KeySort[tau],
    "TauErr"   -> KeySort[tauSE],
    "PlotGamma"-> plotGamma,
    "PlotTau"  -> plotTau,
    "Table"    -> table
  |>
];



MakeMomemtsNormalized[pairs__List,moms_] := 
  Module[{items = pairs, curves = {}, labels = {}},  
   Do[Module[{R = items[[i, 1]], z = items[[i, 2]], x, n, \[Mu]}, 
     x = Select[Flatten@Normal@N@{R}, NumberQ];
     If[Length[x] >= 2, 
     n = Length[x]; \[Mu] = Mean[x];
      AppendTo[curves, 
       Table[{z, m,Mean[(x/\[Mu])^m]},{m,moms}]];
      AppendTo[labels, 
       "z=" <> ToString@NumberForm[z, {8, 3}]<>", <R>=" <> ToString@NumberForm[\[Mu], {8, 3}]]
       ]], {i, Length[items]}];
   If[curves === {}, 
    Return@Graphics[{}, Frame -> True, FrameLabel -> {"m, gamma(m)"}]];
    curves
     ];
   


If[!ValueQ[z0], z0 = 1.];
If[!ValueQ[cycles], cycles = 240];
If[!ValueQ[numdata], numdata = 1000];

MakeZRuns[nz_]:=
Module[{z = 2.,data = {}, RR},
	Do[
	   RR = RunBatch[z, cycles, numdata];
	   AppendTo[data,{RR,z}];
	   z /= 2.
	   ,{nz}
	];
	MakeTailsNormalized@@data
	];


MakeZMomRuns[nz_, moms_]:=
Module[{z = 2.,data = {}, RR, mu, plot},
	Do[
	   RR = RunBatch[z, cycles, numdata];
	   mu = Mean[RR];
	   AppendTo[data,ParallelTable[{z, m,(RR[[k]]/mu)^m,mu},{m,2,moms},{k,Length[RR]}]];
	   z /= 2
	   ,{nz}
	];
	fits =FitMuAndGammas[data];
	<|
    "fits"     -> fits,
    "curves"    -> data
  |>];


MakeZRuns[10]


runs = MakeZMomRuns[20,50];
res = runs["fits"];
curves =res["curves"];


Print["mu-index",res["muindex"]] ;                        (* Around[\[Mu], \[Sigma]_\[Mu]] *)
res["PlotGamma"]   
res["PlotTau"]
res["Table"];

