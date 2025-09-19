
(* ====================== SQG1.wl (Euclidean tetrad seed) ======================
   Patched on 2025-09-18T05:25:29.686775 by Eily: switch to Euclidean tetrad basis.
   - Work entirely in tetrad indices M=1..4 ≡ (t_E, r̂, \[Theta]̂, φ̂), with δ_{AB}.
   - SeedFromZ now returns SU(2) seed P (4×4) in tetrad basis; Abelian row left 0.
   - Optional Urbantke check included (Euclidean triple-product form), local-only.
============================================================================= *)

If[!ValueQ[$SQGTol], $SQGTol = 1.*^-12];

(* Curvature scale ψ(z) with z = r^3 / r0 (M=1) *)
Clear[psiFromZ];
psiFromZ[z_?NumericQ] := 1./(2. z);

(* Pauli (local copy for seed check) *)
If[!ValueQ[pauli],
  pauli = {
    {{0., 1.}, {1., 0.}},
    {{0., -I}, {I, 0.}},
    {{1., 0.}, {0., -1.}}
  };
];
tBasis[A_] := If[A == 1, I IdentityMatrix[2], pauli[[A - 1]]];

(* Build 2×2 from tetrad column M *)
Pcol2x2[P_, M_Integer?Positive] := Sum[tBasis[A] * P[[A, M]], {A, 1, 4}];

(* Euclidean Fhat in tetrad indices (local for the check) *)
Clear[FhatTetradLocal];
FhatTetradLocal[P_] := Module[{F = ConstantArray[0., {3, 4, 4}], PM, PN, comm},
  Do[
    PM = Pcol2x2[P, M]; PN = Pcol2x2[P, N];
    comm = PM.PN - PN.PM;
    Do[
      F[[i, M, N]] = Re[-I Tr[pauli[[i]].comm]];
      F[[i, N, M]] = -F[[i, M, N]],
      {i, 1, 3}
    ],
    {M, 1, 4}, {N, M + 1, 4}
  ];
  Developer`ToPackedArray[F]
];

(* Fast Euclidean Urbantke (tetrad; for seed check only) *)
Clear[UrbantkeMetricTetradFastLocal];
UrbantkeMetricTetradFastLocal[P_, tol_ : $SQGTol] := Module[{F, FF1, FF2, FF3, G, g},
  F = FhatTetradLocal[P];
  FF1 = F[[2]].F[[3]] - F[[3]].F[[2]];
  FF2 = F[[3]].F[[1]] - F[[1]].F[[3]];
  FF3 = F[[1]].F[[2]] - F[[2]].F[[1]];
  G = -(FF1.F[[1]] + FF2.F[[2]] + FF3.F[[3]])/6.;
  G = Chop[(G + Transpose[G])/2, tol];
  (* det-normalize *)
  Module[{d = Det @ N @ G},
    g = If[NumericQ[d] && d =!= 0, N @ Chop[G/Abs[d]^(1/4), tol], IdentityMatrix[4]];
  ];
  <|"densitized" -> G, "g" -> g|>
];

(* Anchored profile n^A on the tiny loop (Euclidean tetrad; unit spatial) *)
Clear[nProfileTetrad];
nProfileTetrad[\[Theta]_] := {0., -Sin[\[Theta]], Cos[\[Theta]], 0.};

(* Euclidean BH seed in tetrad basis (columns M = {t_E, r̂, \[Theta]̂, φ̂}) *)
Clear[SeedFromZ];
SeedFromZ[z_?NumericQ] := Module[{ψ = psiFromZ[z], α, β, P},
  α = 0.5 Sqrt[ψ];  (* along r̂ *)
  β = -Sqrt[ψ];     (* along \[Theta]̂, φ̂ *)
  P = ConstantArray[0., {4, 4}];
  P[[2, 2]] = α;    (* P^1_{r̂} *)
  P[[3, 3]] = β;    (* P^2_{\[Theta]̂} *)
  P[[4, 4]] = β;    (* P^3_{φ̂} *)
  P[[1, 1]] = I/(α β^2);
  Developer`ToPackedArray @ N @ P
];

(* Optional seed check: Urbantke -> δ_AB (Euclidean) *)
Clear[VerifySeedEuclid];
VerifySeedEuclid[z_?NumericQ, tol_ : 1.*^-10] := Module[{P, U, g, δ},
  P = SeedFromZ[z];
  U = UrbantkeMetricTetradFastLocal[P, tol];
  g = U["g"]; δ = IdentityMatrix[4];
  <|"maxAbsError" -> Max[Abs[g - δ]], "g" -> g|>
];

$SQG1Loaded = True;

(* Back-compat alias *)
Clear[nProfile]; nProfile[\[Theta]_] := nProfileTetrad[\[Theta]];
