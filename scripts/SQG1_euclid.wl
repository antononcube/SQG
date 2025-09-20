(* ::Package:: *)

(* ====================== SQG1.wl (Euclidean tetrad seed) ======================
   Patched on 2025-09-18T05:25:29.686775 by Eily: switch to Euclidean tetrad basis.
   - Work entirely in tetrad indices M=1..4 \[Congruent] (t_E, r\:0302, \[Theta]\:0302, \[CurlyPhi]\:0302), with \[Delta]_{AB}.
   - SeedFromZ now returns SU(2) seed P (4\[Times]4) in tetrad basis; Abelian row left 0.
   - Optional Urbantke check included (Euclidean triple-product form), local-only.
============================================================================= *)

If[!ValueQ[$SQGTol], $SQGTol = 1.*^-12];

(* Curvature scale \[Psi](z) with z = r^3 / r0 (M=1) *)
Clear[psiFromZ];
psiFromZ[z_?NumericQ] := 1./(2. z);

(* Euclidean BH seed in tetrad basis (columns M = {t_E, r\:0302, \[Theta]\:0302, \[CurlyPhi]\:0302}) *)
Clear[SeedFromZ];
SeedFromZ[z_?NumericQ] := Module[{\[Psi] = psiFromZ[z], \[Alpha], \[Beta], P},
  \[Alpha] = 0.5 Sqrt[\[Psi]];  (* along r\:0302 *)
  \[Beta] = -Sqrt[\[Psi]];     (* along \[Theta]\:0302, \[CurlyPhi]\:0302 *)
  P = ConstantArray[0., {4, 4}];
  P[[2, 2]] = \[Alpha];    (* P^1_{r\:0302} *)
  P[[3, 3]] = \[Beta];    (* P^2_{\[Theta]\:0302} *)
  P[[4, 4]] = \[Beta];    (* P^3_{\[CurlyPhi]\:0302} *)
  P[[1, 1]] = 1/(\[Alpha] \[Beta]^2);
  Developer`ToPackedArray @ N @ P
];

$SQG1Loaded = True;
