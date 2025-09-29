(* ::Package:: *)

(* ======================= Initialization Helpers ======================= *)

Clear[Ep];
Ep = Table[
    Table[
      Which[
        j <= 3 && k <= 3, LeviCivitaTensor[3][[i, j, k]],
        j == 4 && k <= 3, KroneckerDelta[i, k],
        j <= 3 && k == 4, KroneckerDelta[i, j],
        True, 0
      ],
      {j, 4}, {k, 4}
    ],
    {i, 3}
  ];

Clear[EtaMixed, BarEtaMixed];
EtaMixed = Table[
    Table[
      Which[
        j <= 3 && k <= 3, LeviCivitaTensor[3][[i, j, k]],
        j <= 3 && k == 4, KroneckerDelta[i, j],
        j == 4 && k <= 3, -KroneckerDelta[i, k],
        True, 0
      ],
      {j, 4}, {k, 4}
    ],
    {i, 3}
  ];

BarEtaMixed = Table[
    Table[
      Which[
        j <= 3 && k <= 3, LeviCivitaTensor[3][[i, j, k]],
        j <= 3 && k == 4, -KroneckerDelta[i, j],
        j == 4 && k <= 3, KroneckerDelta[i, k],
        True, 0
      ],
      {j, 4}, {k, 4}
    ],
    {i, 3}
  ];

Aeta := Table[EtaMixed[[a]], {a, 1, 3}];
BAeta := Table[BarEtaMixed[[a]], {a, 1, 3}];

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

(* Sample invocation (uncomment to test interactively) *)

LR = SPQBuildLR[];
z0 = 1.;
sol = SPQSolvePQ[z0, LR, 100, 1.*^-16];
Print["Residual: ", sol["Residual"]];
