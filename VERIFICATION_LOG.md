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
| 5 | spherical harmonics | 🔄 in progress | 5.1 + 5.2 explicit forms started | `check_5_1_2.py` |
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
2. **Chapters 6, 7** — never symbolically verified (OCR + labels only).
3. Chapter 3 is only partially covered (3.2); the rest could be swept.
   (Chapters 0–2 are definitional and not amenable to symbolic checking.)

## Chapter 5 (spherical harmonics) — in progress (2026-08-27)

Reference: `Y_{lm} = mpmath.spherharm(l,m,ϑ,φ)` (VMK = standard Condon–Shortley;
validated in Ch.4). `check_5_1_2.py` verifies (all PASS): 5.2.1
(e^{imφ}√·P_l^m), 5.1.10, 5.1.11 (conjugation), 5.1.6 (orthonormality), the
Rodrigues/differential forms 5.2.2/3/6, hypergeometric 5.2.23, the D-relation
5.2.37, and Gegenbauer 5.2.39.

- **Chapter-wide OCR fix**: `Y_{\text{Im/lm/im}}` (16 spots, incl. section
  titles) → `Y_{lm}`; §5.1.7 `\mathfrak{Z}{lm}`→`\mathfrak{Z}_{lm}`,
  `\vect{3}`→`\mathfrak{Z}`; eq 5.1.9 `\delta_{l}`→`\delta_{ll'}`.
- **Open flags (need scan)**: **5.2.27** off by exactly `|m|!` for |m|≥2 —
  a missing `1/|m|!` (cf. 5.2.23); **5.1.14** `u_lm` differs by √2 (real-harmonic
  normalization — book writes `½(Y+Y*)` but the RHS carries `√((2l+1)/2π)`).
- **§5.2 power-series & hypergeometric cluster** (`check_5_2_series.py`,
  all PASS): the trig-θ/2 series 5.2.9–5.2.14 and hypergeometric
  5.2.24/25/26/31/33/34 verified vs mpmath.spherharm; scan-read the rest
  (printed pp.135–137). Fixes applied: 5.2.9 `∑_*`→`∑_s`; 5.2.17 `(-1)^{(4+m)/2}`
  →`(-1)^{(s+m)/2}` and `(l+ϑ)!!`→`(l+s)!!`; 5.2.20 `∑_ρ`→`∑_s` + √-unwrap of
  `(cosϑ)^l`; 5.2.25 denom `|m||2|m|`→`|m|!2^{|m|}`; 5.2.26 `(\cos^{ϑ/2})`→
  `(\cos ϑ/2)`; 5.2.28 `(\cot^{ϑ/2})`,`(\sin^{ϑ/2})` de-garble; 5.2.33
  `2^{|m||m|!}`→`2^{|m|}|m|!`.
- **Book misprints** (restored + annotated in-source): **5.2.27, 5.2.28**
  omit the `1/|m|!` prefactor (cf. 5.2.23/5.2.25) — added; **5.1.14** `u_lm`
  coefficient `2π`→`4π` (`½(Y+Y*)=Re(Y)` needs 4π, matching v_lm).
- **§5.2 remaining cluster** (`check_5_2b.py`): verified 5.2.21, 5.2.22 (|Y|²),
  5.2.29, 5.2.30, 5.2.32, 5.2.38 vs mpmath.spherharm; fixes applied — 5.2.21
  √-unwrap + `∑_0`→`∑_s` + exponent `l+m/2-s`→`(l+m-s)/2`; 5.2.36 `ρ/2`,`s/2`
  →`ϑ/2`; 5.2.38 `2|m|l!`→`2^{|m|}l!` + √-unwrap; 5.2.40 `1/r`→`1/r^l`; 5.2.18
  `(l+m-1)!`→`(l+m-1)!!`. Chapter now builds clean (fig fixed to `fig5_a`).
- **§5.4 / §5.7 / §5.8** (`check_5_4_7_8.py`, all PASS): symmetry 5.4.1/2/4–6/8;
  recursions 5.7.1–5.7.9; derivative relations 5.8.4/5/6 — verified vs
  mpmath.spherharm. Fixes: **5.4.3** LHS `Y_{lm}`→`Y_{l̄ m}` (l̄=−l−1) + heading;
  **5.7.9** `\sin²\forall`→`\sin²ϑ` and denominator `(2l+1)(2l+3)`→`(2l+1)(2l-3)`
  (confirmed numerically); **5.8.2** garbled `C_{lm+μ}^{lm+μ}`→CG
  `\clebsch{l}{m}{1}{μ}{l}{m+μ}`; **5.8.4** `Y(v,φ)=imY(θ,φ)`→`Y(ϑ,φ)=imY(ϑ,φ)`;
  **5.8.7** `∂/∂ϑ²`→`∂²/∂ϑ²`; **5.8.17** `j_i`→`j_l` (×3) + orphaned `\left.`/`\mid`
  bracket repair; §5.8.1 lead-in `\vect{Y}_{\text{Im}}`→`Y_{lm}`. 5.4.7/5.4.9
  (negative-θ) hold under VMK's `(sinϑ)^m` continuation, not reproducible via
  spherharm's cos-even branch — noted, not flagged.
- **Still deferred in §5.2**: 5.2.15/16 (correct prefactor is
  `√((2l+1)/4π(l+m)!(l-m)!)·l!`, verified — but the book's radical
  organization needs a careful re-read before editing), 5.2.18 (partial fix
  only; my series encoding still off), 5.2.19, 5.2.35/36 (|z|=1 exponential 2F1,
  scan-matched, not numerically reproducible). §5.3 onward untouched.
- **§5.13** (`check_5_13.py` + `check_5_13_1.py`, all PASS): every live form
  verified vs mpmath.spherharm to machine precision.
  - §5.13.2 (eqs 5.13.1–5.13.5, |m|=0..4 via Legendre P_l) and §5.13.3
    (eqs 5.13.6–5.13.11, |m|=l..l−5 via trig) — all correct as printed, no fixes.
  - §5.13.1 (explicit l≤5 table, currently inside a `\begin{comment}` block):
    all 35 entries correct in both the sin^k and multiple-angle forms. OCR label
    garbles fixed: `\sin4\theta`→`\sin4\vartheta` & `e^{-iS\varphi}`→`e^{-i3\varphi}`
    (Y_{4-3}); `\sin3\theta`→`\sin3\vartheta` (Y_{5+5}); `Y_{5+9}`→`Y_{5+3}`;
    `Y_{\delta+1}`/`\sin\phi`→`Y_{5+1}`/`\sin\vartheta`. **Restored dropped
    Y_{5-2}** (OCR omission) as conjugate partner of Y_{5+2}, with in-source NOTE.
- **§5.14** (`check_5_14.py`, all PASS): Y and ∂_ϑY at ϑ=0,π/2,π,±nπ
  (eqs 5.14.1–5.14.8) verified vs mpmath.spherharm (Y at the special ϑ directly;
  ∂_ϑY via high-precision finite differences at the poles / centered at π/2;
  5.14.8 cross-checked through the verified 5.8.5c recurrence to avoid the
  past-the-pole continuation branch). **Fix: 5.14.6** phase
  `(-1)^{l+(m+1)/2}`→`(-1)^{(l+m+1)/2}` — OCR pulled l out of the fraction
  numerator; scan (PDF p.171) confirms the whole l+m+1 is over 2. Silent OCR fix.
- **§5.10** (`check_5_10_1.py` + `check_5_10_2.py`, all PASS): sums involving Y.
  - §5.10.1 (sums over m, eqs 5.10.1–5.10.5) verified to ~1e-30. **Fix 5.10.3**
    lower limit `\sum_{m=-1}^{l}`→`\sum_{m=-l}^{l}` (OCR; scan p.163 confirms −l).
  - §5.10.2 (sums over l, eqs 5.10.6–5.10.16) verified by truncating the |t|<1
    series: Gegenbauer finite sum (6), ₀F₁ product (7), ₂F₁ gen.fn (8), ₂F₁
    product (9), Laguerre gen.fn (10), Rayleigh plane-wave j_l (11/12/13),
    Bessel addition theorem (15, z=y_l), bilinear J_m (16). 5.10.14 is the x=y=t
    boundary case of 5.10.15 (direct sum diverges, terms ~l^{m-1/2}); confirmed
    via the exact identity RHS_14 = RHS_15|_{x=y=t}. OCR fixes (scan p.163-164):
    5.10.7 `o_1(;m+1;.)`→`{}_0F_1`; 5.10.11-14 `j(t)`→`j_l(t)` (spherical Bessel);
    5.10.16 `j i(t)`→`j_l(t)` and split merged exponential
    `e^{it cosϑ₁cosϑ₂ · e^{im(φ₁-φ₂)}}`→`e^{it cosϑ₁cosϑ₂} e^{im(φ₁-φ₂)}`.
- **§5.11** (`check_5_11.py`, all PASS): generating functions. Term-wise
  convergent 5.11.1 (1/R^{2m+1}, both |t|≶1), 5.11.2 (1/R, m=0), 5.11.3
  ([(1+R)²−t²]^{−m}/R) verified by direct summation. 5.11.4/5/7 are
  distributional (Θ step-fn + LHS singularity at ϑ=ψ, not term-wise convergent);
  5.11.7 (m=0 Legendre) confirmed via Abel summation (r→1); 5.11.4/5 scan-matched
  (PDF p.165). OCR fixes: 5.11.1 `t^{+m+1}`→`t^{l+m+1}`; 5.11.2 `1/1^{l+1}`→
  `1/t^{l+1}`; 5.11.3 LHS bracket `\left.\mid…\right]`→`\left[…\right]` and stray
  `\|_{1}`→`l!` (factor pinned numerically); 5.11.7 `\theta`→`\Theta` (×2).
- **§5.9** (`check_5_9.py`, all PASS): integrals involving Y (φ-integral done
  analytically → 1-D θ-quadratures; CG/3j from sympy).
  - §5.9.1 solid-angle 5.9.1–5.9.5 (orthonormality + Gaunt/three-Y). OCR fixes:
    5.9.1 `Y_{1m}`/`\delta_{10}`→`Y_{lm}`/`\delta_{l0}`; 5.9.4 first factor
    `Y_{l₂m₁}`→`Y_{l₁m₁}` and CG J-index `{l₂}{0}`→`{l₃}{0}`; 5.9.5 first factor
    `Y_{l₁m₂}`→`Y_{l₁m₁}` and 3j row `m₁m₃m₃`→`m₁m₂m₃` (scan p.161).
  - §5.9.2 Fourier 5.9.6–5.9.9: distributional/operator (δ(q−k), L̂, ∇×L̂) — not
    numerically verified; notation fixed `j l`/`j i (qr)`→`j_l(qr)` (×3).
  - §5.9.3 θ-integrals 5.9.10–5.9.13. OCR: 5.9.10 `\delta l'`→`\delta_{ll'}`;
    5.9.12 factor `(l+m)!/(l-m)!` moved *inside* the sqrt (scan p.162).
    **BOOK MISPRINT 5.9.13**: closed-form coefficients printed as bare (l₂+m₂),
    (l₁+m₁) — the un-normalized P_l^m form — corrected to
    √((2l+1)(l²−m²)/(2l−1)) required for normalized Y_lm; derived from the
    verified 5.8.6 recurrence, confirmed for general m₁,m₂; annotated in-source.
- **§5.6** (`check_5_6.py`, all PASS): expansions of products of Y (CG from
  sympy). Verified 5.6.9 (Clebsch–Gordan series), 5.6.10 (inverse), 5.6.11
  (three-Y), 5.6.14/5.6.15 (irreducible tensor products), 5.6.17 (iterated
  {Y₁⊗…⊗Y₁} coupling) vs mpmath.spherharm to ~1e-17. §5.6.1 (completeness /
  Parseval / Dirac notation) is definitional. OCR fixes (scan p.157): 5.6.9
  garbled 2nd CG `C_{l₂m₁l₂m a}`→`\clebsch{l₁}{m₁}{l₂}{m₂}{L}{M}`; 5.6.10 `m₃`→
  `m₂`; 5.6.11 `(4x)²`→`(4π)²`; 5.6.14 `\vect{Y}_{l₃}`→`\vect{Y}_{l₂}`; 5.6.15
  stray trailing `,(` removed.
- **§5.3** (`check_5_3.py`, all PASS): integral representations, verified by
  mpmath.quad/quadosc. §5.3.1: Mehler–Dirichlet 5.3.2/5.3.3 (m=0 only — m≥1 is a
  finite-part integral), 5.3.4. §5.3.2: 5.3.5–5.3.10 (both ± signs). §5.3.3:
  improper 5.3.12/5.3.13, 5.3.14 (cosϑ>0), |Y|² Bessel 5.3.15. OCR fixes:
  **5.3.9, 5.3.10** `(1+m)!`→`(l+m)!` (l→1; numerically confirmed); `\item(c)`→
  `\item` (redundant literal label); 5.3.15 `\subsubsection{(c) …}`→`\item …`
  (broke the enumerate). **Domain note (in-source)**: the reciprocal-power forms
  5.3.6/5.3.8/5.3.10 are valid only for cosϑ>0 — for cosϑ<0 the
  ∫dψ/(a+b cosψ) branch flips sign and the RHS returns −Y_lm; the direct forms
  5.3.5/5.3.7/5.3.9 hold for all ϑ. 5.3.1 (m-fold indefinite) and 5.3.11
  (change of variables) left as structural, not numerically checked.

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
