(* ::Package:: *)

(* Wigner (small) d-function d^J_{M M'}(beta) in algebraic form -- VMK convention.

   CONVENTION.  Mathematica's built-in WignerD[{j,m,mp},beta] is the TRANSPOSE of
   the Varshalovich-Moskalev-Khersonskii (VMK) function:

       WignerD[{j, m, mp}, beta] == d^J_{mp,m}(beta)   (VMK)
                                 == (-1)^(m-mp) d^J_{m,mp}(beta)
                                 == d^J_{m,mp}(-beta).

   e.g. WignerD[{1/2,1/2,-1/2},beta] = +Sin[beta/2], whereas VMK has -Sin[beta/2].
   So VMK's d^J_{m,mp}(beta) is obtained by SWAPPING the two projections:

       wignerD[j,m,mp] := WignerD[{j, mp, m}, beta].

   Verified against the explicit Wigner sum formula for j=1/2,1,3/2,2,5/2 and all
   m,m'.  Run  "Get" on this file, then use wignerD / wignerDMatrix below, or run
       WolframKernel -noprompt -script wigner_d.wl
   to see the demo and the self-test result.

   Half-integers must be exact (1/2, 3/2, ...), not machine reals, to stay algebraic. *)

(* VMK-convention small-d, algebraic in the symbol beta (\[Beta]) *)
wignerD[j_, m_, mp_, b_ : \[Beta]] := FullSimplify[WignerD[{j, mp, m}, b]];

(* Full VMK d^j(beta); rows m and cols m' run j, j-1, ..., -j *)
wignerDMatrix[j_, b_ : \[Beta]] := Module[{o = Range[j, -j, -1]},
   FullSimplify[Outer[wignerD[j, #1, #2, b] &, o, o]]];

(* --- independent Wigner explicit sum formula (VMK), used only for the self-test --- *)
refD[j_, m_, mp_, b_] := Sqrt[(j + m)! (j - m)! (j + mp)! (j - mp)!]*
   Sum[(-1)^k Cos[b/2]^(2 j - 2 k + m - mp) Sin[b/2]^(2 k - m + mp)/
       (k! (j + m - k)! (j - mp - k)! (k - m + mp)!),
     {k, Max[0, m - mp], Min[j + m, j - mp]}];

selfTest[] := Module[{Js = {1/2, 1, 3/2, 2, 5/2}, vals},
   vals[j_] := Range[j, -j, -1];
   AllTrue[Js, Function[j, AllTrue[vals[j], Function[m, AllTrue[vals[j],
      Function[mp, FullSimplify[wignerD[j, m, mp, \[Beta]] - refD[j, m, mp, \[Beta]]] === 0]]]]]]];

(* Demo / validation printed when this file is run as a script *)
Print["self-test wignerD == VMK formula (j=1/2..5/2): ", selfTest[]];
Print["d^{1/2}_{1/2,-1/2}(beta) = ", wignerD[1/2, 1/2, -1/2]];
Print["d^{2}_{1,0}(beta)        = ", wignerD[2, 1, 0]];
Print["d^1(beta) = ", wignerDMatrix[1]];
