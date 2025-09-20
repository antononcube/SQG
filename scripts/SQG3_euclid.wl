(* ::Package:: *)

(* ====================== SQG3_euclid.wl ====================== *)
(* Gaussian Fourier forcing sampled along SDE trajectories.       *)
(* ============================================================ *)

$HistoryLength = 0; (* Avoid kernel history bloat. *)

(* import SQG1, SQG2 into this file's namespace, as if they were one file *)
(* Resolve base folder robustly on any kernel (wolframscript or remote) *)
base = Quiet@Check[DirectoryName[$InputFileName], $Failed];
If[!StringQ[base] || base === "", base = Directory[]];
If[!FileExistsQ[FileNameJoin[{base, "SQG1_euclid.wl"}]], base = FileNameJoin[{base, "scripts"}]];

Get[FileNameJoin[{base, "SQG1_euclid.wl"}]];
Get[FileNameJoin[{base, "SQG2_euclid.wl"}]];

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
    chunkSize_: 10^4,
    flushFun_: Automatic
  ] :=
  Module[{\[Theta]s, d\[Theta], P, f\[Theta], sol, k, R, stride = thinStride, chunk = chunkSize,
          buffer = {}, localCollect = {}, flush},

    If[!IntegerQ[stride] || stride <= 0, stride = 1];
    If[!IntegerQ[chunk]  || chunk  <= 0, chunk  = 10^4];

    flush = If[flushFun === Automatic,
      Function[list, localCollect = Join[localCollect, list]],
      flushFun
    ];

    \[Theta]s = Subdivide[0., 2 Pi, nSteps - 1];
    d\[Theta] = \[Theta]s[[2]] - \[Theta]s[[1]];
    P  = SeedFromZ[z];
    R  = RfromP[P, z, s];

    For[k = 1, k <= nSteps, k++,
      f\[Theta]  = fFun[\[Theta]s[[k]]];
      sol = StepSDE20[P, f\[Theta], d\[Theta], stabEvery, k, tol];
      P   = sol["P"];

      If[Mod[nSteps - k, stride] == 0,
        R = RfromP[P, z, s];
        AppendTo[buffer, R];
        If[Length[buffer] >= chunk,
          flush[buffer];
          buffer = {};
        ];
      ];
    ];

    If[buffer =!= {}, flush[buffer]];

    If[flushFun === Automatic, localCollect, Null]
  ];



(* Run directory helper. *)
$RunsDir = "runs";
Clear[NewRunDir];
NewRunDir[] := Module[{files, indices, next},
  If[Not[DirectoryQ[$RunsDir]], CreateDirectory[$RunsDir]];
  files   = FileNames[RegularExpression["\\d+"], $RunsDir];
  indices = ToExpression /@ (FileNameTake[#, -1] & /@ files);
  next    = If[indices === {}, 0, Max[indices] + 1];
  FileNameJoin[{ $RunsDir, IntegerString[next, 10, 4]}]
];

(* Main driver: parallelises over seed list and streams results to disk. *)
Clear[RunNewSimulation];
RunNewSimulation[
    z_?NumericQ,
    nSteps_: 512,
    M_: 6,
    K_: 12,
    sigma_: 0.15,
    s_: 1,
    stabEvery_: 15,
    thinStride_: 10,
    tol_: $SQGTol,
    useeds_: Range[100],
    chunkSize_: 10^4
  ] := Module[{rundir, params, paramsFile, rThinFile, tempFileForSeed, tempFiles,
    destHandle},

  rundir = NewRunDir[];
  CreateDirectory[rundir];
  Print["Run: ", rundir];

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
    "useeds" -> useeds
  |>;
  paramsFile = FileNameJoin[{rundir, "params.txt"}];
  Export[paramsFile, params];
  Export[FileNameJoin[{rundir, "date_begin.txt"}], Now];

  rThinFile = FileNameJoin[{rundir, "rthin.bin"}];
  If[FileExistsQ[rThinFile], DeleteFile[rThinFile]];

  tempFileForSeed = Function[seed,
    FileNameJoin[{rundir, "rthin_" <> IntegerString[seed, 10, 4] <> ".bin"}]
  ];

  DistributeDefinitions[MakeUProfile, SampleWTrajectory, chunkSize, tempFileForSeed];

  ParallelMap[
    Function[seed,
      Module[{f, tempFile, writer},
        tempFile = tempFileForSeed[seed];
        If[FileExistsQ[tempFile], DeleteFile[tempFile]];
        f = MakeUProfile[M, K, sigma, seed];
        writer = Function[list,
          If[list === {}, Return[]];
          Module[{fh = OpenAppend[tempFile, BinaryFormat -> True]},
            BinaryWrite[fh, Developer`ToPackedArray @ N[list], "Real64"];
            Close[fh];
          ];
        ];
        SampleWTrajectory[z, f, s, nSteps, stabEvery, tol, thinStride,
          chunkSize, writer];
        Print["Completed SampleWTrajectory for seed: ", seed, "."];
        Null
      ]
    ],
    useeds,
    ProgressReporting -> True
  ];

  destHandle = OpenWrite[rThinFile, BinaryFormat -> True];
  tempFiles = Table[tempFileForSeed[seed], {seed, useeds}];
  Scan[
    Function[file,
      If[FileExistsQ[file],
        Module[{fh = OpenRead[file, BinaryFormat -> True]},
          BinaryWrite[destHandle, BinaryReadList[fh, "Real64"], "Real64"];
          Close[fh];
        ];
        DeleteFile[file];
      ]
    ],
    tempFiles
  ];
  Close[destHandle];

  Export[FileNameJoin[{rundir, "date_end.txt"}], Now];
];


(* Demo run parameters (adjust as needed). *)
On[Assert];
RunNewSimulation[
  1.,      (* z *)
  1024,    (* nSteps *)
  6,       (* M *)
  12,      (* K *)
  1.1,     (* sigma *)
  1,       (* s *)
  1,       (* stabEvery *)
  10,      (* thinStride *)
  $SQGTol, (* tol *)
  Range[1] (*seeds*)
];
