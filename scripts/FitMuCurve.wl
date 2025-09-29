Clear[FitMuCurve];
(* finite numeric *)
isFiniteQ[expr_] := NumericQ[expr] && FreeQ[expr, _DirectedInfinity | Indeterminate];

Options[FitMuCurve] = {
  "Degree" -> 2,           (* polynomial degree in x = -Log[z] *)
  "CompressByZ" -> True,   (* average repeats at same z and weight by n/var *)
  "SnapZ" -> False,        (* snap z (or x) to bins to merge float noise *)
  "SnapDigits" -> 8        (* digits for snapping if SnapZ->True *)
};

FitMuCurve[quads_, OptionsPattern[]] := Module[
  {deg = OptionValue["Degree"], compress = OptionValue["CompressByZ"],
   snap = OptionValue["SnapZ"], digs = OptionValue["SnapDigits"],
   x, flat, data, pairs, byX, binned, X, Y, W, lm, coeffs, se, r2, aic, bic,
   xMin, xMax, fitFun, cov, secondFun, secondSEFun, xGrid, y2Band, plot, tbl},

  (* 1) flatten and extract {x=-log z, y=log mu}, keep finite *)
  flat = Flatten[quads, 1];
  data = Reap[
    Do[
      With[{xx = -Log[flat[[i, 1]]], yy = Log[flat[[i, 4]]]},
        If[ isFiniteQ[xx] && isFiniteQ[yy], Sow[{xx, yy}] ]
      ],
      {i, Length[flat]}
    ]
  ][[2, 1]];
  If[data === {} , Return[$Failed, Module]];

  (* optional: snap x to integer grid of 'digs' digits to merge float noise *)
  If[snap === True, data[[All, 1]] = Round[data[[All, 1]], 10.^-digs]];

  (* 2) compress by repeated x if requested *)
  If[compress,
    byX = GroupBy[data, First -> Last];
    binned = KeyValueMap[
      Function[{xx, ys},
        Module[{n = Length[ys], ybar = Mean[ys], s2 = If[Length[ys] > 1, Variance[ys], 0.]},
          {xx, ybar, If[s2 > 0., n/s2, n/10.^6], n, Sqrt[s2/Max[1, n - 1]]}  (* {x, ȳ, w, n, stderr} *)
        ]
      ],
      byX
    ] // Values // SortBy[#, First]&;
    X = binned[[All, 1]]; Y = binned[[All, 2]]; W = binned[[All, 3]];
    ,
    X = data[[All, 1]]; Y = data[[All, 2]]; W = ConstantArray[1., Length[Y]];
    binned = Transpose[{X, Y, W, ConstantArray[1, Length[Y]], ConstantArray[0., Length[Y]]}];
  ];

  If[Length[Y] < (deg + 1) || Variance[X] == 0., Return[$Failed, Module]];

  (* 3) weighted polynomial fit y = Sum_{k=0..deg} a_k x^k *)
  lm = Quiet @ LinearModelFit[Transpose[{X, Y}], Table[x^k, {k, 0, deg}], x, Weights -> W];

  coeffs = lm["BestFitParameters"];        (* a0..adeg *)
  se     = lm["ParameterErrors"];
  r2     = lm["RSquared"];  aic = lm["AIC"];  bic = lm["BIC"];

  fitFun = lm[#] &;
  xMin = Min[X]; xMax = Max[X];

  (* 4) second derivative y''(x) and its 1-σ from covariance *)
  cov = lm["CovarianceMatrix"];            (* (deg+1)x(deg+1) *)
  secondFun[x_] := Sum[k (k - 1) coeffs[[k + 1]] x^(k - 2), {k, 2, deg}];
  (* grad wrt parameters: g_j(x) = d y'' / d a_j *)
  With[{d = deg},
    secondSEFun[x_] := Module[{g},
      g = Table[If[j >= 3, (j - 1)(j - 2) x^(j - 3), 0.], {j, 1, d + 1}];
      Sqrt[g . cov . g]
    ]
  ];

  (* 5) build plot with error bars (if compressed) and fit curve *)
  xGrid = Subdivide[xMin, xMax, 400];
  plot = Show[
    If[compress,
      ListPlot[
        Transpose @ { binned[[All, 1]], Around @@@ Transpose[{binned[[All, 2]], binned[[All, 5]]}]},
        PlotMarkers -> Automatic, PlotStyle -> Black, PlotTheme -> "Scientific"
      ],
      ListPlot[Transpose @ {X, Y}, PlotMarkers -> Automatic, PlotStyle -> Black]
    ],
    Plot[fitFun[x], {x, xMin, xMax}, PlotStyle -> {Red, Thick}],
    Frame -> True, FrameLabel -> {"x = -log z", "y = log μ"}, ImageSize -> 520
  ];

  (* 6) present coefficient table *)
  tbl = Grid[
    Prepend[
      Table[{Subscript[a, k], coeffs[[k + 1]], se[[k + 1]]}, {k, 0, deg}],
      {"coef", "value", "stderr"}
    ],
    Frame -> All
  ];

  <|
    "DataBinned"      -> binned,          (* {x, mean, weight, n, stderr} *)
    "Model"           -> lm,
    "Coefficients"    -> coeffs,
    "CoeffSE"         -> se,
    "R2"              -> r2, "AIC" -> aic, "BIC" -> bic,
    "FitPlot"         -> plot,
    "Second"          -> secondFun,
    "SecondSE"        -> secondSEFun,
    "SecondAt"        -> (x \[Function] Around[secondFun[x], secondSEFun[x]]),
    "XRange"          -> {xMin, xMax}
  |>
];
