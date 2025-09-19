(* load suite *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "scripts", "SQG1_euclid.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "scripts", "SQG2_euclid.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "scripts", "SQG3_euclid.wl"}]];

(* simple seed check *)
res = VerifySeedEuclid[10.];
err = res["maxAbsError"] /. _Missing -> 1.;

If[NumericQ[err] && err < 1.*^-8,
  Print["PASS: VerifySeedEuclid[10.] maxAbsError = ", err]; Exit[0],
  Print["FAIL: VerifySeedEuclid[10.] maxAbsError = ", err]; Exit[1]
];
