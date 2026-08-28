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
| 4 | Wigner D-functions | ✅ complete | **§4.3–4.19 verified** (2026-08-28); no open flags | `check_4_3.py` … `check_4_16_17.py`, `check_4_18.py`, `check_4_19.py` |
| 5 | spherical harmonics | ✅ complete | **all machine-checkable eqs in §5.1–5.17 verified** (2026-08-27/28); no open flags | `check_5_1_2.py`, `check_5_2_series.py`, `check_5_2b.py`, `check_5_3.py`, `check_5_4_7_8.py`, `check_5_5.py`, `check_5_6.py`, `check_5_9.py`, `check_5_10_1.py`, `check_5_10_2.py`, `check_5_11.py`, `check_5_12.py`, `check_5_13.py`, `check_5_13_1.py`, `check_5_14.py`, `check_5_15.py`, `check_5_16.py`, `check_5_17a.py`, `check_5_17b.py` |
| 6 | spin functions / density matrix | ✅ complete | **§6.1–6.3 verified** (2026-08-28); 5 book misprints fixed | `check_6_1.py`, `check_6_2.py`, `check_6_3.py` |
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

**Chapter 4 complete** (2026-08-26/28). §4.3–4.19 all verified; no open flags.
Next candidates: Chapters 6, 7 (never symbolically verified) or the rest of
Chapter 3.
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
- Sec 4.18 (asymptotics / infinitesimal rotations) — `check_4_18.py`: a fast
  mpmath `dJ` (self-tested vs `wigner_d` to 1e-16, dps=200 for large-J
  cancellation) verifies 4.18.1 (large-J cos asymptotic, confirmed O(J^-3/2)
  by error-ratio 8.5 ≈ 4^1.5 over a J×4 step), 4.18.2 (Bessel limit),
  4.18.3/4 (β→0 and
  π−β→0 expansions, confirmed O(β^4) residual), and 4.18.5–8 (infinitesimal
  rotations = −i⟨JM|J·n|JM'⟩, exact). **Three book misprints fixed** (all
  scan-confirmed, printed p.116, annotated in-source):
  - 4.18.2: Bessel order `J_{M-M'}(Jβ)` → `J_{M'-M}(Jβ)` (differ by (−1)^{M−M'}).
  - 4.18.3: leading `(β/2)^μ` → `(sin β/2)^μ` (the exact threshold factor; only
    then is the printed O(β²) bracket accurate to O(β⁴)).
  - 4.18.4: leading `((π−β)/2)^ν` → `(sin((π−β)/2))^ν` (same, near β=π).
- Sec 4.19 / Table 4.2 (other authors' conventions) — `check_4_19.py`:
  reference data, mostly convention-dependent. Verified the definitional row
  (authors with e^{+i} rotation operators ⇒ D = D_book(−α,−β,−γ)) and D_book
  unitarity. OCR fix: **`e^{∓iγ Ĵ_z}` had been OCR'd as `e^{∓i† Ĵ_z}`
  (`\dagger`)** in 4 rotation-operator cells (Edmonds/Dolginov/Rose[30]/
  Gel'fand). Phase/transpose rows (Gel'fand −i vs +i, Yutsis, Berestetskii)
  depend on each author's own handedness/rotation-sense and are left as
  published (documented in the script). (Aside: cited name "Brink and
  Satcher" should read "Satchler" — not touched, outside equation scope.)

1. **Chapter 4 continued** — 4.4 (symmetries), 4.6–4.7 (sums / addition
   theorems), 4.8 (recursions), 4.10–4.12 (orthogonality / integrals),
   4.14–4.15 (characters). All machine-checkable via `wigner_d`.
2. **Chapter 7** — never symbolically verified (OCR + labels only). (Chapter 6
   completed 2026-08-28.)
3. Chapter 3 is only partially covered (3.2); the rest could be swept.
   (Chapters 0–2 are definitional and not amenable to symbolic checking.)

## Chapter 6 (spin functions) — complete (2026-08-28)

All three sections verified as finite-dimensional linear algebra; the physics is
correct, most defects were OCR. Spin-S / T_LM(S) operators built from Sec 2.4/2.6.

- **§6.1 (arbitrary spin)** — `check_6_1.py`: single-particle 6.1.5, 6.1.11–15,
  6.1.18, 6.1.22, 6.1.27, 6.1.30–35, 6.1.49–51 for S=½,1,3⁄2; two-particle
  6.1.53–65 + 6.1.34 for (S₁,S₂)=(½,½),(1,½),(1,1). Uses a **Racah-formula 6j**
  (sympy `wigner_6j` raises spuriously on some valid half-integer 6j, e.g.
  {1 1 1;1 ½ ½}). **Book misprints:** 6.1.60 & 6.1.65 6j `{S₁ S₂ L; S S S}` is
  identically 0 for mixed spins (non-integer triad sum) → Q_L eigenvalue would
  vanish; correct `{S₁ S₂ S; S₂ S₁ L}` (eigenvalue-proven). 6.1.58 phase
  `(-1)^{2S+L₂}` → `(-1)^{S₁+S₂+S′+L₂}` inside the sum (Edmonds 7.1.7; the
  particle-1 analog 6.1.57 with `(-1)^{S₁+S₂+S+L₁}` is correct). 6.1.46 book
  `ρ²±ρ` → `ρ²=ρ`. OCR: pervasive lowercase s→S, merged S_χ subscripts, garbled
  CG in 6.1.50, e_x′→e_z′, 6.1.58 6j lower-row S′/L₂ swap.
- **§6.2 (S=½)** — `check_6_2.py`: basis/helicity spinors, product expansions
  6.2.9/10, rotated bases 6.2.17/18 vs D^{½}/U^{½}, helicity 6.2.22–31,
  spin-direction 6.2.39–41, polarization 6.2.50–52, contravariant transform
  6.2.46. ~46 OCR fixes (χ scanned as x/X, ½ as ⅓/¼, S_z as S_x, γ→α in 6.2.17
  exponents, φ→ψ in 6.2.22/23, dropped conjugate star in 6.2.46 — scan-checked).
- **§6.3 (S=1)** — `check_6_3.py`: T_{2M} via Eq 2.6.4, Q_{ik} via 2.6.7;
  product expansions 6.3.20–23 (middle + cartesian), action tables 6.3.24–33
  (emitted and reconciled), rotated/helicity forms 6.3.36–52, ⟨T_{2M}⟩=…Y_{2M}
  (6.3.62). ~37 OCR fixes, dominated by χ_x/χ_z and Q_xx/Q_zz/Q_zx diagonal-index
  garbles and S_x↔S_z operator labels.

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
- **§5.12** (`check_5_12.py`, all PASS): asymptotics, verified as limits
  (relative error vanishes at the stated order). 5.12.1 (leading O(1/l)), 5.12.3
  (bound), 5.12.4/5.12.5 (small-ϑ / near-π, O(δ⁴)), 5.12.6/5.12.7 (near-π/2,
  O(δ⁴)), 5.12.8 (McDonald Bessel), 5.12.9 (J_m(lϑ), rel err O(1/l)).
  **BOOK MISPRINT 5.12.2**: printed second-term sign "+" makes the "more exact"
  formula worse than the leading 5.12.1 for m≥1 (error stays O(1/l)); with "−"
  it is genuinely O(1/l²), 2–27× better (improvement grows with m). Corrected +
  annotated. **OCR 5.12.6**: phase `(-1)^{(l+m)/2}`→`(-1)^{(l±m)/2}` (scan p.166
  shows the ±; the + form fails for Y_{l,−m} with odd m).
- **§5.5** (`check_5_5.py`, all PASS): behaviour under coordinate transforms.
  Angle relations 5.5.2 (inversion), 5.5.6/7/8 (π-rotations x/y/z), 5.5.9 (z by
  χ), 5.5.11 (equatorial refl.), 5.5.12 (meridian refl.); general rotation 5.5.1
  via Wigner D (passive R⁻¹); infinitesimal 5.5.10 (d/dω vs ladder bracket);
  translation 5.5.3. OCR fixes: 5.5.1 `D^{\prime}`→`D^l` and `Y(ϑ,ψ)`→`Y(ϑ,φ)`;
  5.5.4 `d²`→`a²`; 5.5.9 `e^{-imx}`→`e^{-imχ}`. **BOOK MISPRINT 5.5.3**: the
  coefficient printed as `[4π(2l+1)(2l'-2l+1)/(2l'+1)]^{1/2}` is wrong for l'≥2
  (coincides with the correct value only at l'=1); the regular-solid-harmonic
  addition theorem gives `[4π(2l'+1)!/((2l+1)!(2l'-2l+1)!)]^{1/2}` — factorials
  dropped and fraction inverted. Corrected + annotated (verified to 1e-15).
- **§5.16** (`check_5_16.py`, all PASS): bipolar & tripolar harmonics (CG/3j/9j
  from sympy). Bipolar: 5.16.4 (Σ|Bip|²), 5.16.6 (inversion phase), 5.16.7/5.16.8
  (product CG series — one 9j), 5.16.9/5.16.10 (L=0 scalar product). Tripolar:
  5.16.14 (Σ|Trip|²), 5.16.17/5.16.18 (product CG series — TWO 9j, verified to
  3e-19), 5.16.19 (L=0 tripolar scalar, 3j). Definitional: 5.16.1/2/3/5/11-13/15/16.
  OCR fixes (scan p.174-175): **5.16.11** garbled CG `C_{l₂m₂m₂m₃m₃}^{l₂₃}`→
  `C^{LM}_{l₁m₁ l₂₃m₂₃} C^{l₂₃m₂₃}_{l₂m₂ l₃m₃}`; **5.16.7/5.16.8** B-subscript
  `l₁'l₂'L_1''l₂''L''`/missing l₂''→`l₁'l₂'L'l₁''l₂''L''`; **5.16.17/5.16.18**
  B superscript missing l₂ and subscript `l'`/`L_1''` garbles → full index list.
  The B-coefficient *formulas* (9j structure) were correct as printed.
- **§5.17** (`check_5_17a.py` + `check_5_17b.py`, all PASS): expansions of
  two-vector functions. Scalar group (17 identities, reduced to 1-D Legendre via
  the addition theorem (Y_l·Y_l)=(2l+1)/4π P_l(cosω₁₂)): 5.17.9, (r₁·r₂)ⁿ
  5.17.10–13, exp 5.17.14, power-series 5.17.16/17, Green's fns 5.17.19/21/22,
  r/1/r³/1/r⁵ 5.17.23–25, and every hypergeometric aₗⁿ coefficient 5.17.26–33 —
  all correct as printed to ~1e-30. Bipolar group: 5.17.34 (spherical wave, z=n),
  5.17.35 (rᴸY_LM), 5.17.36 (Y_LM/r^{L+1}), 5.17.37 (r_μ/r³), 5.17.38/39 (rᴺY_LM).
  **Note (5.17.39)**: the generic-N coefficient degenerates at Γ-pole lattice
  points — (L-N)/2 ≤0 integer (polynomial cases, → 5.17.35) and (L+N+3)/2 ≤0
  integer (very-negative N); exact wherever both Γ's are regular (10 (N,L) pairs
  verified). OCR fixes: 5.17.7 `l₃`→`l₂`; 5.17.14/5.17.19 `j l`/`j(x)`→`j_l(x)`;
  5.17.34 `\dot{l}_1`→`j_{l₁}`. 5.17.1-8 (setup) / 5.17.18 (δ-fn) definitional.
- **§5.15** (`check_5_15.py`, all PASS): zeros of Y_lm and ∂_ϑY_lm. All
  cos²ϑ_α (5.15.2, m=l..l−5) and cos²ϑ_β (5.15.4, m=l..l−4) formulas confirmed by
  plugging the predicted cos²ϑ back and checking Y (resp. ∂_ϑY) vanishes there
  (~1e-30 / 1e-15); interior zero counts (5.15.1) verified = l−|m|. **OCR 5.15.4**:
  the m=±(l-3) row was mislabeled `m=±(l-4)` (duplicate of the next row) and
  carried a spurious `; 0` — scan p.173 shows `m=±(l-3) … ; 1(l≥5)`; the checker
  verified that formula as l−3. Also `\partial\theta`→`\partial\vartheta` in the
  count text. 5.15.3 (large-l approx θ_α) and Table 5.1 left as reference.

- **§5.2 stragglers resolved** (`check_5_2b.py`, all PASS; scan pp.148,150-151):
  the 6 previously-flagged forms are now verified to ~1e-15/1e-40.
  - **5.2.15/16**: OCR mis-scoped the radical (only (l+m)!(l-m)! under it; l! and
    (cos/sin th/2)^{2l} outside) AND the book carries a **spurious |m|!** — correct
    prefactor l! sqrt((l+m)!(l-m)!) (matches Wigner d^l_{m,0}). Both corrected.
  - **5.2.18**: even-branch denominator `s!(sin th)^s`→`s!!(sin th)^s` (OCR; scan
    confirms double factorial). **5.2.19**: radical mis-scoped (wraps only the
    (l-m)!/(l+m)! ratio; (sin th)^m and sum outside) + `L+m-s`→`l+m-s`; m≥0 form.
  - **5.2.35/36**: **BOOK MISPRINT** in leading sign — printed −i/π yields exactly
    −Y_lm; correct is +i/π. 5.2.36 clinches it (its 2F1 arg has |z|=1/(2 sin th)<1
    near π/2 — convergent series, no branch ambiguity — yet still −Y with −i).
    Corrected to +i/π (verified 6e-41).

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

## Lessons learned (methodology — from the Chapter 5 pass, 2026-08-27/28)

Distilled from verifying all of §5.1–5.17 and cracking the six §5.2 stragglers.
These are transferable to the remaining chapters.

**1. Verification-first, always.** Build the numeric checker for a section
*before* editing the `.tex`. It is the only reliable way to tell an OCR slip
(fix silently) from a genuine book misprint (fix *and* annotate). Guessing from
the scan alone repeatedly mislabeled both.

**2. Push past the first non-trivial order.** Several book errors are *exactly
right at the lowest order* and only diverge higher up — they survive casual
checks precisely because dipole-level (l=1) sanity tests pass:
  - 5.5.3 translation coefficient: exact at l'=0,1, wrong for l'≥2.
  - 5.10.14, 5.11.3: right at low order, structural error above.
Verify at l up to ~6–8 and m across the full range, not just a token case.

**3. The discrepancy's *shape* names the bug.** This was the single most useful
diagnostic. Compute ratio (or residual) vs the exact Y over several (l,m,θ):
  - **Constant ±1, independent of l,m** → global sign / branch / domain issue,
    NOT a coefficient typo. (5.2.35/36 sign; 5.3.6/8/10 domain; the −Y tell.)
  - **θ-dependent ratio** → structural error inside the summand/coefficient.
    (5.2.18 s!→s!! showed up this way.)
  - **Right at low l, wrong higher** → wrong *functional form* of a coefficient;
    fit it from a small linear system, then identify. (5.5.3.)
  - **Selective failure** (only Y_{l,−m}, only odd m, only one L-parity) →
    a sign-of-m / parity tracking error. (5.12.6 (l±m)/2; 5.15.4 mislabel.)
  - **A phase that must be real but isn't** (½-integer exponent) → a mis-split
    fraction. (5.14.6 (l+m+1)/2.)

**4. Know when a formula is *supposed* to break.** Not every failure is an
error. Recognize and test the *regular* regime instead of chasing the tail:
  - Generic-N / hypergeometric closed forms have Γ-pole lattice points
    (5.17.39: (L−N)/2 or (L+N+3)/2 a non-positive integer).
  - Distributional / Abel-summable series (5.11.4/5/7, 5.10.14) — verify by
    Abel summation (r→1⁻) or as the limit of a neighbouring convergent case.
  - Domain restrictions (5.3.6/8/10 need cosθ>0; the reciprocal-power branch
    flips sign otherwise).
  - Analytic-continuation subtleties can *look* unresolvable — find the
    sub-case where they evaporate. (5.2.36's |z|=1/(2sinθ)<1 near π/2 is a plain
    convergent series and proved the sign error that 5.2.35's |z|=1 obscured.)

**5. Hard cases are usually *stacked* errors.** The stragglers resisted because
each combined an OCR slip *and* a book misprint (5.2.15/16 = mis-scoped radical
+ spurious |m|!). Fix one at a time and re-measure; a partial fix leaves a
partial residual that looks like a different bug.

**6. Exploit identities to collapse the check.** Reducing dimensionality makes
high-precision verification cheap and robust:
  - Addition theorem (Y_l·Y_l)=(2l+1)/4π P_l(cosω₁₂) turns every two-vector
    *scalar* expansion (§5.17) into a 1-D Legendre identity.
  - Doing the φ-integral analytically turns solid-angle integrals (§5.9) into
    1-D θ-quadratures.
  - Derive a doubted coefficient from an *already-verified* neighbour (5.9.13's
    K_l came straight from the verified 5.8.6 recurrence).

**7. Cross-check corrections against theory, not just numerics.** When a fix is
substantial, confirm it also drops out of a known representation — 5.2.15/16's
corrected l!√((l+m)!(l-m)!) prefactor matches the Wigner d^l_{m,0} form; 5.5.3's
corrected [4π(2l'+1)!/((2l+1)!(2l'−2l+1)!)]^{1/2} is the solid-harmonic addition
theorem. Two independent routes to the same answer ≫ one numeric fit.

**8. Recurring OCR signatures in this corpus** (grep for these first):
  - `\forall`/`\exists`/`\theta` in math where `\vartheta` is meant (ϑ↔∀).
  - single `!` for `!!` (double factorial); `s!`↔`s!!`.
  - `j l`, `j i`, `\dot{l}` → `j_l` (spherical Bessel subscript lost).
  - `\sqrt{…}` wrapping too much (a trailing trig factor or l! swept under the
    vinculum) — the most common §5.2 signature.
  - `(4x)`→`(4π)`; capital `L`→`l`; duplicated case labels (m=±(l−4) twice);
    garbled CG / recoupling-symbol indices (subscripts dropped or swapped).

**9. Tooling / workflow.**
  - `mpmath.spherharm(l,m,θ,φ)` = VMK Y (Condon-Shortley); `sympy.physics.wigner`
    for CG/3j/6j/9j (slow — cache results, it is the bottleneck; drop dps and
    truncation for infinite sums, r1/r2≈0.3 converges by l≈40).
  - Scan: PDF page = printed page + 13; `pdftoppm -r 200-300` then crop with PIL.
  - Build hygiene: delete `ChapN.aux` before a `\includeonly{ChapN}` rebuild;
    note the `sed`/`perl` includeonly rewrite can silently no-op (whole book
    builds, 375pp) — 0 errors still validates the chapter.
  - Run slow checkers in the background (`-u` for progressive output), wait with
    Monitor; never `git checkout -- VMK.tex` (clobbers uncommitted figure work);
    stage only the specific files for each section's commit.

_Last updated: 2026-08-28 (Chapter 5 complete)._
