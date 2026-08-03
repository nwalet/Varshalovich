(* ::Package:: *)

(* U-function U^J_{M M'}(omega; Theta, Phi) in algebraic form -- VMK convention.

   U^J_{M M'}(omega; Theta, Phi) = <J M| Exp[-I omega n.J] |J M'> is the matrix of
   a rotation by the angle omega about the axis
       n(Theta, Phi) = (Sin[Theta] Cos[Phi], Sin[Theta] Sin[Phi], Cos[Theta]),
   in the convention of Varshalovich-Moskalev-Khersonskii (VMK), Sec. 4.5.
   It reduces to the small Wigner d-function for a rotation about the y-axis:
       U^J_{M M'}(omega; Pi/2, Pi/2) == d^J_{M M'}(omega).

   The closed form below has a REMOVABLE singularity at Theta = 0 (there v = 0
   and the 1/v^2 terms cancel analytically): use the symbolic result, or keep
   Theta != 0 for purely numerical evaluation.  Pass half-integer J, M, M' as
   exact rationals (1/2, 3/2, ...) so the result stays algebraic.

   Verified against the matrix exponential Exp[-I omega n.J] for j = 1/2 .. 2
   (see selfTest[] at the bottom; run this file as a script to print it). *)

Ufunction[j_, m_, mp_, w_ : \[Omega], th_ : \[CapitalTheta], ph_ : \[CapitalPhi]] :=
  Block[{uu, v = Sin[w/2] Sin[th], ampp = Abs[m + mp], maxs},
   uu = Cos[w/2] - If[m + mp >= 0, I, -I] Sin[w/2] Cos[th];
   maxs = If[m + mp >= 0, Min[j - m, j - mp], Min[j + m, j + mp]];
   FullSimplify[
    (-I v)^(2 j - ampp) uu^ampp Exp[-I (m - mp) ph]*
     Sqrt[(j + m)! (j - m)! (j + mp)! (j - mp)!]*
     Sum[If[m + mp >= 0,
            1/(s! (s + ampp)! (j - m - s)! (j - mp - s)!),
            1/(s! (s + ampp)! (j + m - s)! (j + mp - s)!)] (1 - v^-2)^s,
       {s, 0, maxs}]]];

(* Full U^J(omega; Theta, Phi); rows M and cols M' run J, J-1, ..., -J *)
UfunctionMatrix[j_, w_ : \[Omega], th_ : \[CapitalTheta], ph_ : \[CapitalPhi]] :=
  Module[{o = Range[j, -j, -1]},
   Table[Ufunction[j, mm, mmp, w, th, ph], {mm, o}, {mmp, o}]];

(* --- independent reference for the self-test: <J M| Exp[-I w n.J] |J M'> --- *)
Jz[j_] := DiagonalMatrix[Table[m, {m, j, -j, -1}]];
Jp[j_] := Module[{ms = Table[m, {m, j, -j, -1}], n}, n = Length[ms];
   Table[If[i + 1 == k, Sqrt[j (j + 1) - ms[[k]] (ms[[k]] + 1)], 0], {i, n}, {k, n}]];
Jx[j_] := (Jp[j] + Transpose[Jp[j]])/2;
Jy[j_] := (Jp[j] - Transpose[Jp[j]])/(2 I);
refU[j_, m_, mp_, w_, th_, ph_] :=
  Module[{nn = {Sin[th] Cos[ph], Sin[th] Sin[ph], Cos[th]}, U},
   U = MatrixExp[-I w (nn[[1]] Jx[j] + nn[[2]] Jy[j] + nn[[3]] Jz[j])];
   U[[j - m + 1, j - mp + 1]]];

selfTest[] := Module[
   {pts = {{0.7, 1.1, 0.3}, {2.0, 0.4, 1.7}, {1.3, 2.5, -0.9}}, worst = 0, vals},
   vals[j_] := Range[j, -j, -1];
   Do[Do[Do[Do[
       worst = Max[worst, Abs[N[Ufunction[j, m, mp, p[[1]], p[[2]], p[[3]]] -
             refU[j, m, mp, p[[1]], p[[2]], p[[3]]]]]],
      {p, pts}], {mp, vals[j]}], {m, vals[j]}], {j, {1/2, 1, 3/2, 2}}];
   worst];

(* Demo / validation printed when this file is run as a script *)
Print["self-test  max|Ufunction - Exp[-I w n.J]|  (j=1/2..2) = ", selfTest[]];
Print["U^{1/2}_{1/2,1/2}(w; Th, Ph) = ", Ufunction[1/2, 1/2, 1/2]];
Print["U(w; Pi/2, Pi/2) for j=1  (== d^1(w)) = ",
  FullSimplify[UfunctionMatrix[1, \[Omega], Pi/2, Pi/2]]];
