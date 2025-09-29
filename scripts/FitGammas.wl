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

  (* Group by m → {{x,y},...} *)
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
  tau  = AssociationThread[validM, (m \[Function] (m*muindex + gamma[m])) /@ validM];
  tauSE = AssociationThread[validM, (m \[Function] Sqrt[(m*muSE)^2 + (gammaSE[m])^2]) /@ validM];

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
      {"m", "\[Gamma](m)", "StdErr(γ)", "\[Tau](m)", "StdErr(τ)", "R^2(γ-fit)"}
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
