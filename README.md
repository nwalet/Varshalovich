# Quantum Theory of Angular Momentum — electronic edition

A re-typeset, fully searchable LaTeX edition of

> D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii,
> **Quantum Theory of Angular Momentum**,
> World Scientific, Singapore, 1988.

The original is the standard reference on the subject, but it is available only
as a scan: not searchable, uneven in image quality, and awkward to read on
anything but paper. This project rebuilds it as real LaTeX — selectable text,
live cross-references, a proper index, and vector diagrams — and adds machinery
to **verify the formulas** rather than trusting the OCR.

The original was made Open Access in 2021 through the sponsorship of SCOAP³ and
is licensed **CC BY 4.0**, which permits this adaptation. See
[License and attribution](#license-and-attribution).

---

## What is here

| | |
|---|---|
| **14 chapters** | `Chap0.tex` (introduction) through `Chap14.tex`, plus `references.tex`, `glossary.tex` |
| **327 diagrams** | redrawn as standalone TikZ in `tikz_files/`, built to `images/` |
| **83 Python programs** | 73 verification checkers, 3 table generators and assorted helpers in `scripts/` |
| **2 Jupyter notebooks** | interactive generation of algebraic and numerical table entries |
| **5 LaTeX packages** | `diagrams.sty` and friends, for the graphical-method notation |

### Chapters

| | | | |
|---|---|---|---|
| 1 | Elements of Vector and Tensor Theory | 8 | Clebsch–Gordan Coefficients and 3*jm* Symbols |
| 2 | Angular Momentum Operators | 9 | 6*j* Symbols and the Racah Coefficients |
| 3 | Irreducible Tensors | 10 | 9*j* and 12*j* Symbols |
| 4 | Wigner *D*-Functions | 11 | The Graphical Method |
| 5 | Spherical Harmonics | 12 | Sums Involving Recoupling Coefficients |
| 6 | Spin Functions | 13 | Matrix Elements of Irreducible Tensor Operators |
| 7 | Tensor Spherical Harmonics | 14 | Computer Generation of Numerical and Algebraic Expressions |

Chapter 14 is **new to this edition**. It documents the computer-algebra
methods below; the other thirteen follow the original.

---

## Verification, not just transcription

Re-typesetting a 500-page book of formulas by OCR introduces errors. A wrong
index or a dropped factor is invisible on the page but fatal in use, so most
formulas here are checked mechanically against independent implementations.

`scripts/` holds 73 checkers, one per section or family of identities. Each
evaluates the book's formula and an independent reference at many numerical
points and reports any disagreement:

```bash
python3 scripts/check_9_11.py       # 6j algebraic tables vs sympy
python3 scripts/check_10_13.py      # the 12j identities against each other
```

This has found genuine defects — both OCR slips and, occasionally, misprints
carried over from the printed book. Where a formula legitimately differs from
the naive reading, the checker records why.

### Generating the tables

The algebraic tables of §8.12, §9.11 and §10.11 are not transcribed but
**generated**, and can be regenerated for any entry — including values beyond
the ranges the book prints:

```bash
python3 scripts/gen_8_12_cg_tables.py  --b 1 --beta=0 --k=0
python3 scripts/gen_9_11_6j_tables.py  --d 1 --m=0 --n=0
python3 scripts/gen_10_11_9j_tables.py --alpha 3/2 --beta 3/2 --gamma 0 \
                                       --lam=1/2 --mu=1/2 --nu=0
```

Each prints the entry in the book's own notation and as a SymPy expression.
Omit the row and column indices to dump a whole table.

The method is described in Chapter 14. In outline: each table fixes one small
angular momentum, which makes the relevant Racah sum finite; a change to
integer-valued variables ($p=a+\alpha$, $q=a-\alpha$ for the Clebsch–Gordan
coefficients, the triangle deficits $u=s-2a$ etc. for the 6*j* symbols) turns
every factorial argument into an explicit integer offset, so the factorials
collapse to rising and falling factorials and the algebra stays rational — no
general-purpose simplification anywhere. The 9*j* symbols reduce to a short sum
over three 6*j* symbols and reuse the same engine.

Every generated entry has been checked against SymPy at both the level of the
sum and the level of the regrouping into printed form:

| Table | Cells | Test points | Mismatches |
|---|---|---|---|
| §8.12 Clebsch–Gordan, *b* ≤ 5 | 505 | 31 185 | 0 |
| §9.11 6*j*, *d* ≤ 4 | 284 | 61 344 | 0 |
| §10.11 9*j*, all 11 tables | 438 | 4 344 | 0 |

### Notebooks

Two notebooks give the same functions interactively, with menus where
`ipywidgets` is installed and plain edit-and-run cells where it is not:

- **`vmk_algebraic_tables.ipynb`** — algebraic entries for the Clebsch–Gordan,
  6*j* and 9*j* tables, as typeset formulas, SymPy expressions or LaTeX source.
- **`vmk_numeric_tables.ipynb`** — exact numerical values as a signed root of a
  rational fraction plus a decimal, in the form the book's numerical tables
  use, covering the Clebsch–Gordan coefficients and the 6*j*, 9*j* and both
  kinds of 12*j* symbol. A coefficient that vanishes reports *which* condition
  failed rather than silently returning zero.

---

## Building the book

Requires a TeX Live installation with **XeLaTeX** (the fonts and
`ucharclasses` make this mandatory — pdfLaTeX will refuse), plus Python 3 with
SymPy for the scripts.

```bash
# 1. build the figures (once; ~330 standalone TikZ files)
cd tikz_files && make -j8 && cd ..

# 2. build the book
xelatex VMK
bibtex  VMK                                   # biblatex, bibtex backend
makeindex VMK.idx                             # subject index
makeindex -s vmksym.ist -o VMK.snd VMK.sdx    # symbol index (-s is required)
xelatex VMK
xelatex VMK
```

The `-s vmksym.ist` flag is not optional: without it the symbol index still
typesets but loses its Latin/Greek/other group headings.

To rebuild only some chapters, uncomment and edit the `\includeonly` line in
`VMK.tex`.

### Repository layout

```
VMK.tex            master document (preamble, includes)
preamble.tex       title page, license, preface
Chap0–14.tex       the chapters
diagrams.sty       macros for the Chapter 11–12 graphical method
tikz_files/        327 standalone TikZ figures + Makefile
images/            figure PDFs (generated; not tracked)
scripts/           73 verification checkers, 3 table generators, helpers
figures/           Mathematica/Asymptote sources for a few 3-D figures
orig/              the original scanned PDF, for comparison
```

---

## Status

The text of all fourteen chapters is complete and indexed. Verification
coverage is extensive but not uniform: the vector-addition and recoupling
coefficients (Chapters 8–10) are checked exhaustively, while other chapters are
checked section by section. §10.13 in particular has 17 of its 36 labelled
equations under test.

Corrections and additional checkers are welcome. If you find a discrepancy with
the printed original, please open an issue quoting the equation number — and
note whether the printed book agrees with the formula or with the checker,
since a few known misprints originate in the 1988 text.

---

## License and attribution

This edition is licensed under the
[Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/)
(CC BY 4.0) — full text in [`LICENSE`](LICENSE), attribution and statement of
changes in [`NOTICE`](NOTICE). The license permits use, sharing, adaptation, distribution and
reproduction in any medium or format, provided you give appropriate credit,
link to the license, and indicate if changes were made.

It is based on the work originally published by World Scientific Publishing Co.
Pte. Ltd., Singapore, 1988 — ISBN 978-9971-5-0107-5 (hardcover),
978-9971-5-0996-5 (paperback), 978-981-4415-49-1 and 978-981-4578-28-8
(ebook). Original copyright © 1988 World Scientific Publishing Co. Pte. Ltd.
The ebook was converted to Open Access in 2021 through the sponsorship of
**SCOAP³** (Sponsoring Consortium for Open Access Publishing in Particle
Physics), under the same CC BY 4.0 license.

**Changes made in this edition:** the text has been re-typeset from the scan in
LaTeX; all diagrams have been redrawn as vector graphics; the index has been
rebuilt and split into subject and symbol indexes; formulas have been checked
mechanically and OCR errors corrected; and Chapter 14, on the computer
generation of the tables, has been added. This edition is not endorsed by the
original authors or publisher.

Re-typeset by **Niels R. Walet**. AI tools were used to assist with text
recognition, formatting and verification.
