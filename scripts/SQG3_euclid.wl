(* ::Package:: *)

(* ::Input:: *)
(* ======================SQG3 . nb====================== *)
(* Neutrality uses n_\[Mu] P^{0\[Mu]} (P stored with upper \[Mu]) . Gaussian Fourier u(\[Theta]) held fixed during solve and sampling . *)
(* ====================================================*)


(* The following two lines are only needed if launching in a notebook environment. *)

SetDirectory@NotebookDirectory[]; (* Go to the directory this file is in. *)
$HistoryLength=0; (* Don't keep Mathematica history. *)


(* import SQG1, SQG2 into this file's namespace, as if they were one file *)
Get["SQG1.wl"]
Get["SQG2.wl"]

If[!ValueQ[$SQGTol],$SQGTol=1.*^-12];
If[!ValueQ[$SQGZ0],$SQGZ0=1.0];

(*Gaussian Fourier profile*)
Clear[MakeUProfile]
MakeUProfile[K_:12,sigma_:0.15,seed_:Automatic] := BlockRandom[
  Module[{a, b, xi},
    a = RandomVariate[NormalDistribution[0,sigma],K];
    b = RandomVariate[NormalDistribution[0,sigma],K];
    xi = RandomVariate[NormalDistribution[0,sigma]];
    Function[\[Theta],xi +Sum[a[[k]] Sin[k \[Theta]]+b[[k]] Cos[k \[Theta]],{k,1,K}]]
  ];
  RandomSeeding -> seed
];

SampleWTrajectory[z_?NumericQ,uFun_,s_:+1,nSteps_:512,stabEvery_:16,tol_:$SQGTol,thinStride_,rThinAppendFun_] := Module[{\[Theta]s,d\[Theta],P,u\[Theta],sol,k},
\[Theta]s=Subdivide[0.,2 Pi,nSteps-1];d\[Theta]=\[Theta]s[[2]]-\[Theta]s[[1]];
P=SeedFromZ[$SQGZ0];
For[k=1,k<=nSteps,k++,u\[Theta]=uFun[\[Theta]s[[k]]];
sol=StepSDE20[P,u\[Theta],d\[Theta],stabEvery,k,tol];
P=sol["P"];
If[Mod[nSteps - k, thinStride] == 0, rThinAppendFun[CosmoInvariants[P,z,s]["R"]]]
];
];

(* RunOneRealization[z_?NumericQ,nSteps_:512,K_:12,sigma_:0.15,s_:+1,seed_:1,stabEvery_:16,tol_:$SQGTol]:=Module[{u=MakeUProfile[K,sigma,seed],W},
W=SampleWTrajectory[z,u,s,nSteps,stabEvery,tol];
<|"xi"->xi,"W"->W|>]; *)

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
      K_ : 12,
      sigma_ : 0.15,
      s_ : 1,
      stabEvery_ : 15,
      thinStride_ : 10,
      tol_ : $SQGTol,
      useeds : Range[100]
    ] := Module[{rundir, params, paramsFile, rThinFile, rThinFileHandle, rThinAppendFun},

  rundir = NewRunDir[];
  Print["Run:", rundir];

  paramsFile = FileNameJoin[rundir, "params.txt"];
  params = <|
    "z" -> z,
    "nSteps" -> nSteps,
    "K" -> K,
    "sigma" -> sigma,
    "s" -> s,
    "stabEvery" -> stabEvery,
    "tol" -> tol,
    "thinStride" -> thinStride,
    "useeds" -> useeds
  |>;
  Export[paramsFile, params];
  Print["Params:", params];

  Export[FileNameJoin[rundir, "date_begin.txt"], Now];

  rThinFile = FileNameJoin[rundir, "rthin.bin"];
  rThinFileHandle = OpenWrite[rThinFile, BinaryFormat->True];
  rThinAppendFun = Function[r, WithLock[rThinFileHandle, BinaryWrite[rThinFileHandle, N[r], "Real64"]]];

  (* let's do multiple u! *)
  ParallelMap[
    Function[seed,
      Module[{u},
        u = MakeUProfile[K, sigma, seed];
        SampleWTrajectory[z,u,s,nSteps,stabEvery,tol,thinStride,rThinAppendFun];
      ]
    ],
    useeds,
    ProgressReporting -> True
  ];

  Close[rThinFileHandle];

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

RunNewSimulation[1.,512, 12, 0.15, 1, 15, 10,$SQGTol,Range[100]]







(* ====================== Euclidean patch (SQG3) ======================
   Patched on 2025-09-18T05:25:29.686775 by Eily
   - No Urbantke in the dynamics; default metric is Euclidean δ.
   - This file intentionally remains tetrad-agnostic; generators work in Euclidean space.
====================================================================== *)
If[!ValueQ[$EuclidMetric], $EuclidMetric = IdentityMatrix[4]];
