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
| 4 | Wigner D-functions | 🔄 in progress | **4.3, 4.4 done**; 4.5–4.15 remain | `check_4_3.py`, `check_4_4.py` |
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

**Now in progress: Chapter 4** (2026-08-26).
- Sec 4.3 (explicit d-forms) — `check_4_3.py`, all 18 forms 4.3.2–4.3.23
  verified; four source errors fixed (see flags).
- Sec 4.4 (symmetries of d/D) — `check_4_4.py`: d-relations 4.4.1, the full
  31-relation D array 4.4.2, periodicity 4.4.4–4.4.5, and special cases
  4.4.6–4.4.8 all verified; the badly OCR-garbled 4.4.2 array was rebuilt
  from the scan (see flags).

## Gaps / next up

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
