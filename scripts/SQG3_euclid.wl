(* ::Package:: *)

(* ::Input:: *)
(* ======================SQG3 . nb====================== *)
(* Gaussian Fourier u(\[Theta]) held fixed during solve and sampling . *)
(* ====================================================*)

$HistoryLength=0; (* Don't keep Mathematica history. *)

(* import SQG1, SQG2 into this file's namespace, as if they were one file *)
Get[FileNameJoin[DirectoryName[$InputFileName], "SQG1_euclid.wl"]];
Get[FileNameJoin[DirectoryName[$InputFileName], "SQG2_euclid.wl"]];

If[!ValueQ[$SQGTol],$SQGTol=1.*^-12];
If[!ValueQ[$SQGZ0],$SQGZ0=1.0];

(*Gaussian Fourier profile*)
Clear[MakeUProfile];
MakeUProfile[M_:6,K_:12,sigma_:0.15,seed_:Automatic] := BlockRandom[
  Module[{a, b, xi, m},
    a = RandomVariate[NormalDistribution[0,sigma],{M,K}];
    b = RandomVariate[NormalDistribution[0,sigma],{M,K}];
    xi = RandomVariate[NormalDistribution[0,sigma],M];
    Function[\[Theta],
    Table[xi[[m]] +Sum[a[[m,k]] Sin[k \[Theta]]+b[[m,k]] Cos[k \[Theta]],{k,1,K}],{m,1,M}]]
  ],
  RandomSeeding -> seed
];

SampleWTrajectory[z_?NumericQ,fFun_,s_:+1,nSteps_:512,stabEvery_:16,tol_:$SQGTol,thinStride_,rThinAppendFun_] := 
Module[{\[Theta]s,d\[Theta],P,f\[Theta],sol,k,R},
  \[Theta]s=Subdivide[0.,2 Pi,nSteps-1];d\[Theta]=\[Theta]s[[2]]-\[Theta]s[[1]];
  P=SeedFromZ[z];
  R=RfromP[P,z,s];
  Print["(init) z=",z,", R=",R];
  For[k=1,k<=nSteps,k++,
    f\[Theta]= fFun[\[Theta]s[[k]]];
    sol=StepSDE20[P,f\[Theta],d\[Theta],stabEvery,k,tol];
    P=sol["P"];
    (* Quit[]; *)
    If[Mod[nSteps - k, thinStride] == 0, 
    R = RfromP[P,z,s];
    Print["k=",k,", R=",R];
    rThinAppendFun[R]
    ]
  ];
];


(* NewRunDir, no input, outputs a directory *)
$RunsDir = "runs";
Clear[NewRunDir];
NewRunDir[] := Module[{lock, files, next, path},
  (* make runs dir if does not exist *)
  If[Not[DirectoryQ[$RunsDir]], CreateDirectory[$RunsDir]];
  (* make lock file if doesn't exist *)
  lock = FileNameJoin[$RunsDir, ".lock"];
  If[Not[FileQ[lock]], CreateFile[lock]];
  (* make next run dir safely *)
  WithLock[File[lock],
    files = FileNames[RegularExpression["\\d+"], $RunsDir];
    next = Max[Join[{-1}, Map[Composition[ToExpression, FileNameTake], files]]] + 1;
    path = FileNameJoin[$RunsDir, IntegerString[next, 10, 4]];
    CreateDirectory[path]
  ];
  (* return path *)
  path
];

Clear[RunNewSimulation];
RunNewSimulation[
      z_?NumericQ,
      nSteps_ : 512,
      M_ : 6,
      K_ : 12,
      sigma_ : 0.15,
      s_ : 1,
      stabEvery_ : 15,
      thinStride_ : 10,
      tol_ : $SQGTol,
      useeds_ : Range[100]
    ] := Module[{rundir, params, paramsFile, rThinFile, rThinLockFile, rThinAppendFun},

  rundir = NewRunDir[];
  Print["Run: ", rundir];

  paramsFile = FileNameJoin[rundir, "params.txt"];
  params = <|
    "z" -> z,
    "nSteps" -> nSteps,
    "M" -> M,
    "K" -> K,
    "sigma" -> sigma,
    "s" -> s,
    "stabEvery" -> stabEvery,
    "tol" -> tol,
    "thinStride" -> thinStride,
    "useeds" -> useeds
  |>;
  Export[paramsFile, params];

  Export[FileNameJoin[rundir, "date_begin.txt"], Now];

  SetSharedVariable[rThinLock];
  rThinFile = FileNameJoin[rundir, "rthin.bin"];
  rThinFileHandle = OpenWrite[rThinFile, BinaryFormat->True]; (* should this be shared? *)
  rThinAppendFun = Function[r, BinaryWrite[rThinFileHandle, N[r], "Real64"]];

  (* let's do multiple u! *)
  (* ParallelMap[
    Function[seed,
      Module[{f},
        f = MakeUProfile[M, K, sigma, seed];
        SampleWTrajectory[z,f,s,nSteps,stabEvery,tol,thinStride,rThinAppendFun];
        Print["Completed SampleWTrajectory for seed: ", seed, "."];
      ]
    ],
    useeds,
    ProgressReporting -> True
  ];
  Close[rThinFileHandle];
   *)

  
  (* replaced ParallelMap[...] with a serial Do loop *)
    Do[
     Module[{f},
        f = MakeUProfile[M, K, sigma, seed];
        SampleWTrajectory[z,f,s,nSteps,stabEvery,tol,thinStride,rThinAppendFun];
        Print["Completed SampleWTrajectory for seed: ", seed, "."];
      ],
    {seed, useeds}
    ];


  Export[FileNameJoin[rundir, "date_end.txt"], Now];
];

(* Reading data: BinaryReadList["path/to/seed.bin", "Real64"] *)
(* Reading params: Import["path/to/params.txt"] *)

(* todo:
  - Use OptionsPattern. https://reference.wolfram.com/language/tutorial/Patterns.html#14984
  - Fix errors it throws in calculation.
  - Fix error in writing R.
  - Better overview of progress.
*)

On[Assert];
RunNewSimulation[1000., 512, 6, 12, 0.1, 1, 15, 10,$SQGTol, {1}]
