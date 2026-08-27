# Equation Verification Log

Progress record for the **symbolic/numeric verification pass** over VMK
(*Quantum Theory of Angular Momentum*, Varshalovich–Moskalev–Khersonskii).
This is distinct from the OCR-conversion pass tracked in
[CONVERSION_RULES.md](CONVERSION_RULES.md): here each machine-checkable
equation is re-derived and evaluated to confirm the LaTeX matches a correct
identity, and genuine source misprints are flagged.

## Method

- Checkers live in [`scripts/`](scripts/), one file per chapter/section
  (`check_<chap>_<sec>.py`, `check_chap13_*.py`, plus the Chapter-8 topical
  scripts `symmetries_8_4.py`, `special_values_8_5.py`, `recursions_8_6.py`,
  `sums_8_7.py`). Each script's docstring lists the exact equations it covers.
- **Numeric** checks use `sympy.physics.wigner` (`clebsch_gordan`,
  `wigner_6j`, `wigner_9j`, `racah`) evaluated over ranges of j/m values.
- **Symbolic** checks use `sympy` for algebraic identities (R-symbol algebra,
  tensor-product coefficients, operator relations) with free symbols.
- Run one with e.g. `python3 scripts/check_9_1.py`; a clean run prints each
  equation with a PASS/OK marker.
- Convention: when a check fails and the fix is ours (OCR), the `.tex` is
  corrected and the commit says "verify/fix"; when the book itself is wrong,
  the equation is annotated in the `.tex` as a source misprint/omission and
  left mathematically as the *correct* form.

## Status by chapter

| Ch | Topic | Verified? | Sections covered | Scripts |
|----|-------|-----------|------------------|---------|
| 0–2 | prelim / rotations / D-functions | ➖ n/a | not amenable to symbolic checking (definitional) | — |
| 3 | irreducible tensors | ✅ partial | 3.2 tensor-product coefficients | `check_chap3.py` |
| 4 | Wigner D-functions | 🔄 in progress | **4.3–4.17 done**; 4.18–4.19 remain | `check_4_3.py` … `check_4_16_17.py` |
| 5 | — | ❌ **gap** | OCR + headings/labels only | — |
| 6 | — | ❌ **gap** | OCR + headings/labels only | — |
| 7 | — | ❌ **gap** | OCR + headings/labels only | — |
| 8 | CG / 3jm | ✅ extensive | 8.1–8.2 tables, 8.4 symmetry, 8.5 special values, 8.6 recursions (incl. 8.6.1, 8.6.8 Regge), 8.7 sum rules (8.7.1–8.7.6), 8.8 generating fns, 8.10 zero selection rules; accidental 3j zeros (J≤17) | `symmetries_8_4`, `special_values_8_5`, `recursions_8_6`, `sums_8_7`, `check_8_8`, `check_8_10`, `check_cg_tables`, `check_3j_zeros` |
| 9 | 6j / Racah / R-symbol | ✅ complete | 9.1–9.8 | `check_9_1,2,4,5,6,8` |
| 10 | 9j / 12j | ✅ extensive | 10.1–10.2, 10.4–10.9, 10.12–10.13 (incl. 12j(I) & 12j(II)) | `check_10_1,2,4,5,8,9,13` |
| 11 | graphical method | ➖ n/a (diagram chapter) | bra-ket bracket audit only (see below) | — |
| 12 | graphical sums | ➖ n/a (diagram chapter) | — | — |
| 13 | tensor-operator matrix elements | ✅ complete | 13.1; 13.2.1–13.2.9 (nabla, J, L, S, Y_L, scalar/vector products) | `check_chap13_J,L,S,Y,Y2,Y3,Y4,n,nabla,scalar` |

## Where it stopped

Previously: Chapter 13, §13.2.9 (`be21346`, ~2026-08-13). The pass had run
3 → 8 → 9 → 10 → 13.

**Now in progress: Chapter 4** (2026-08-26/27). §4.3–4.17 done; next §4.18
(asymptotics) / §4.19 (other authors' conventions).
- Sec 4.3 (explicit d-forms) — `check_4_3.py`, all 18 forms 4.3.2–4.3.23
  verified; four source errors fixed (see flags).
- Sec 4.4 (symmetries of d/D) — `check_4_4.py`: d-relations 4.4.1, the full
  31-relation D array 4.4.2, periodicity 4.4.4–4.4.5, and special cases
  4.4.6–4.4.8 all verified; the badly OCR-garbled 4.4.2 array was rebuilt
  from the scan (see flags).
- Sec 4.5 (U-matrix in ω,Θ,Φ) — `check_4_5.py`: the explicit forms 4.5.3/4/6,
  the property relations 4.5.17–4.5.22, and special cases 4.5.28–4.5.32
  verified against the validated `u_function` helper; seven errors fixed
  (5 OCR + 2 book misprints, see flags).
- Sec 4.6 (sums of D) — `check_4_6.py`: CG series 4.6.1, projected sums
  4.6.3/5/6, k-fold products 4.6.10, Cayley-Klein form 4.6.13 verified; 3 OCR
  fixes (4.6.5 `M_3`→`M_2`; 4.6.6 spurious `M1M2M` sum + `D_{M_N N}`→`D_{MN}`
  + missing label; 4.6.9 `D_{ƒ(a)M_k N_k}`→`D_{M_k N_k}`).
- Sec 4.7 (addition of rotations) — `check_4_7.py`: unitarity 4.7.4, d-additions
  4.7.7–4.7.10, identical-rotation 4.7.15, and the character sums 4.7.17/4.7.19
  with their angle formulas 4.7.18/4.7.20 all verified; 2 OCR fixes (4.7.16
  `ϖ`→`φ`; 4.7.19 `R_2^{-r}`→`R_2^{-1}`).
- Sec 4.8 (D^J recursions) — `check_4_8.py`: all 21 relations 4.8.1–4.8.21
  (D^{J±1}, D^{J±1/2}, and the same-level M/M'±1 ladders) verified; many OCR
  fixes — the D^{J±1/2} phases `e^{i(α+α)/2}`→`e^{i(α+γ)/2}` and `α±x`→`α±γ`,
  missing `M` subscripts and `D^{J(...)}` garbles in 4.8.4/5/7/8, the 4.8.7 last
  coefficient, and the badly-garbled 4.8.18–4.8.21 (rebuilt from the scan: `B`
  diagonal denominators, `M'`→`M` in 4.8.20/21, and two `sqrt(x/x)=1` coeffs);
  added missing labels 4.8.16/18/19/20.
- Sec 4.9 (differential relations) — `check_4_9.py`: all of 4.9.1–4.9.9 verified;
  1 LaTeX-bracket fix in 4.9.1 and a missing label on 4.9.2.
- Secs 4.10–4.11 (orthogonality/integrals) — `check_4_11.py`: the integral
  identities 4.11.1–4.11.8 verified by β-quadrature (the α,γ integrals are exact
  Kronecker deltas); OCR fixes to 4.11.4 (`D^{J3}_{M2M3'}`→`M3M3'`, CG bottom
  `-M1'`→`-M3'`) and 4.11.5 (`M5'`→`M3'`). 4.10 is prose/definitional.
- Sec 4.12 (invariant summation) — 4.12.2/3/4 verified numerically (character
  integrals → δ and {J1J2J3}). **4.12.5–4.12.9** (multi-R integrals →
  6j / 6j² / three-6j / 9j / 9j²) now fully reconstructed character-by-character
  against the scan (printed pp.97–98, 200 dpi) and confirmed by a clean build.
  Corrections: 4.12.5 spurious `^2` removed (the 6j is NOT squared); 4.12.6 R2
  first D gained its conjugate `*` and the 6j is `{J1J2J3/J1'J2'J3'}^2`
  (all-primed bottom, squared); 4.12.7 R1 2nd D subscript `N2 N3'`→`N2 N2'`;
  4.12.9 R2 2nd D `J6`→`J4` and R4 2nd D `N2 M2`→`N3 M3`. 4.12.8 (the 9j)
  matched as-is. These quadruple-R identities are impractical to evaluate
  numerically; verification is against the printed equations.
- Sec 4.13 (generating functions for d^J) — not machine-verified (would need
  the ξ,μ,ν helpers), but 4 clear OCR garbles fixed against the scan
  (printed p.99): 4.13.6 `-t^{\sin^2}β/2`→`-t\sin^2 β/2` and exponent
  `s+(μ±ν)/2`→`s+(μ+ν)/2`; 4.13.7 exponent `(μ+ν)/2`→`s+(μ+ν)/2`; 4.13.8
  Pochhammer `(μ+ν+1-λ),(λ)_s`→`(μ+ν+1-λ)_s(λ)_s` and exponent
  `+(μ+ν)/2`→`s+(μ+ν)/2`.
- Sec 4.14 (characters χ^J(ω)) — `check_4_14.py`: ~45 identities verified
  (explicit forms 4.14.3–14, properties/periodicity 4.14.16–23, diff. eq
  4.14.26, diff. relation 4.14.30, algebraic 4.14.31–37, orthogonality/
  integrals 4.14.38/40–44, finite sums 4.14.45–51, infinite/generating series
  4.14.52/55–62, particular ω 4.14.63–66, special cases 4.14.67–72) — all
  PASS. Two source errors fixed: 4.14.22 (OCR) and 4.14.62 (book misprint,
  see flags). 4.14.9 is valid only on ω∈[0,π] (noted in the checker).

- Sec 4.15 (generalized characters χ_λ^J) — `check_4_15.py`: ~28 identities
  verified against the definition 4.15.1 (λ-th derivative of χ^J in cos(ω/2),
  cross-checked vs the CG series 4.15.2) — trig series 4.15.2–4, differential
  form 4.15.5, Gegenbauer 4.15.6, Jacobi 4.15.7/8, hypergeometric 4.15.9/10/13/
  15/16, integral rep 4.15.18, symmetries 4.15.19/20, particular ω 4.15.21,
  recursions 4.15.22/23, asymptotics 4.15.25, ODE 4.15.26, orthogonality
  4.15.28 — all PASS. Source fixes: 4.15.21 denominator `(2J−λ)!`→`(2J−λ)!!`
  (OCR); LaTeX de-garbles in 4.15.12 (`(\sin^{ω/2})`→`(\sin ω/2)^λ`), 4.15.13
  (`(\cos^{…})` mangle), 4.15.15/16 (`!11`,`!1`,single-`!` → `!!`). Forms
  4.15.11/12/14 (2F1 with |z|>1 or complex z and negative-integer c) hold only
  under analytic continuation — scan-verified (pp.107–108), not numerically;
  documented in-source.

- Secs 4.16–4.17 (D for particular arguments / particular M,M') —
  `check_4_16_17.py`: 4.16.1–4.16.5 (D at β=0, ±2nπ, ±(2n+1)π, π/2 binomial
  sum) and 4.17.1/2/4/7/8 all verified against D=e^{-iMα}d e^{-iM'γ} (wigner_d)
  and VMK-convention Y (validated via 4.17.1). Phase garbles fixed: **§4.17
  `e^{±i(α±α)/2}`→`e^{±i(α±γ)/2}`** — 4.17.4's four phases numerically
  confirmed (α±γ), 4.17.3's two taken from the scan (p.114, `e^{±i(α−γ)/2}`;
  its general form uses half-integer-degree Y and isn't cleanly evaluable, but
  its M'=±1/2 special case 4.17.4 is). Also 4.17.4 line-1 bracket fix: the
  `cos(β/2)/(J+½)` factor had been pulled into the exponent.

1. **Chapter 4 continued** — 4.4 (symmetries), 4.6–4.7 (sums / addition
   theorems), 4.8 (recursions), 4.10–4.12 (orthogonality / integrals),
   4.14–4.15 (characters). All machine-checkable via `wigner_d`.
2. **Chapters 5, 6, 7** — never symbolically verified (OCR + labels only).
3. Chapter 3 is only partially covered (3.2); the rest could be swept.
   (Chapters 0–2 are definitional and not amenable to symbolic checking.)

## Open flags (book misprints / omissions, annotated in-source)

- **eqs 4.3.3, 4.3.9, 4.3.10** — OCR errors fixed (2026-08-26): 4.3.3 phase
  `(-1)^{J+M'}`→`(-1)^{J+M}`; 4.3.9 inner bracket `(J∓M')`→`(J∓M)`; 4.3.10
  sqrt denominator `(J+M')!`→`(J+M)!` and inner bracket `(J±M')`→`(J±M)`.
  Verified against the scan (printed pp.76–77).
- **eq 4.3.23** — source misprint corrected (2026-08-26): hypergeometric
  argument `-1/cos^2(β/2)`→`+1/cos^2(β/2)`; annotated in-source.
- **eq 4.4.1** — OCR fixes (2026-08-26): row-1 `d_{M'M'}`→`d_{M'M}` and
  `d_{-M'-M'}`→`d_{-M'-M}` (verified vs scan, printed p.79).
- **eq 4.4.2** (D symmetry array) — rebuilt from the scan (2026-08-26): the
  `.tex` array had `M'M'`/`-M'-M'` throughout (should be `M'M`/`-M'-M`), a
  4-argument `D(γ,β,-β,-γ)` garble, a spurious `D^{JJ'}`, and misaligned
  columns. All 31 relations of the corrected array verified numerically.
- **Sec 4.5 fixes** (2026-08-26): eq 4.5.10 dropped prime `J_z`→`J_z'` [OCR];
  eq 4.5.19 truncated arg `(ω;Θ`→`(ω;Θ,Φ)` [OCR]; eq 4.5.20 stray trailing
  `m` [OCR]; eq 4.5.21 `U_{-M'-M'}`→`U_{-M'-M}` [OCR]; eq 4.5.25 garble
  `sin²Θ₁ ..`→`sin²Θ = 1` [OCR]; **eq 4.5.22 book misprint** RHS `π+Φ`→`π-Φ`;
  **eq 4.5.28 book misprint** x-axis `D(π/2,ω,−π/2)=(−i)^{M−M'}d`→
  `D(−π/2,ω,π/2)=(i)^{M−M'}d` (the printed RHS is the inverse rotation).
  Both misprints annotated in-source.
- **eq 4.14.22** — OCR fix (2026-08-26): the transformation J̃=−J−1 was lost;
  `.tex` read `χ^J(R)=−χ^J(R)` (i.e. χ=0). Restored the tilde:
  `χ^{J̃}(R)=−χ^J(R)`, matching the scan (printed p.101). Verified: χ^{−J−1}=−χ^J.
- **eq 4.14.62** — book misprint (2026-08-26; printed p.105). The printed
  denominator `1+t²−4t cos(ω/2)cos(ω′/2)+2t²(cosω+cosω′)` is numerically wrong;
  the correct form (matching the series to 1e-14) restores the (1+t²) factors:
  `(1+t²)²−4t(1+t²)cos(ω/2)cos(ω′/2)+2t²(cosω+cosω′)`. Annotated in-source.
- **eq 9.6.3** — annotated as a source misprint (`e73fc49`).
- **eq 10.8.26** — annotated as a source omission (`9212d0b`).
- **eq 10.9.8** — flagged, then corrected as a source-author fix and verified
  (`ed5745a` → `800f8b5`).
- **eq 10.8.2, 10.8.23, 10.13.34, 10.13.31/32, 10.13.4** — source fixes made
  and verified during the Ch10 pass.
- Chapter 11 (separate audit, 2026-08-26): reduced-matrix-element bracket
  convention confirmed against the scan — round for the standard rank-k
  reduced m.e. (Table 11.14), angle for the generalized ones (Tables 11.9,
  11.15). Book itself carries a `j'`/`m'` misprint in the four-line
  generalized-W-E definition (Chap11.tex ~L1592), reproduced as-is.

_Last updated: 2026-08-26._
