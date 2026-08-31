# Items for review

Issues found while adding index entries to Chapters 0–13. None of these are
indexing problems; they are pre-existing defects that the indexing pass walked
past. Nothing in this list has been changed unless the entry says so.

State verified against a **full 13-chapter build** (`\includeonly` disabled),
two `xelatex` passes, on 2026-08-31. All thirteen chapters are now indexed
(1279 entries: 1092 subject, 194 symbol), including the Introduction.

Current build status: **2 errors, 5 undefined references, 0 undefined citations.**
Both errors and all five references are listed below; none are indexing-related.

---

## A. Build errors (2)

Both are missing figure files in Chapter 11.

```
! LaTeX Error: File `fig11_3_100c' not found.
! LaTeX Error: File `fig11_3_103c' not found.
```

`images/fig11_3_100.pdf` and `images/fig11_3_103.pdf` exist; the `c` variants do
not. Since `fig11_3_90c`, `91c`, `92c` and `96c` all exist, the `c` suffix looks
like the corrected-redraw convention, so these are two redraws not yet done.

**Deliberately not "fixed":** dropping the `c` would compile, but would silently
substitute the *uncorrected* figures. Needs the redraws, or a decision to fall
back.

---

## B. Undefined references (5)

All five render as `??` in the PDF.

| Ref | Used at | Target status | Likely intent |
|---|---|---|---|
| `2:5:36` | [Chap1.tex:1155](Chap1.tex#L1155) | `eq:2:5:36` **exists** in Chap2.tex | `\eqref{eq:2:5:36}` — missing `eq:` prefix |
| `4` | [Chap6.tex:280](Chap6.tex#L280) | `eq:4:1:2` **exists** in Chap4.tex | `\eqref{eq:4:1:2}` — see below |
| `chap8:tab:10` | [Chap8.tex:129](Chap8.tex#L129) | label **is** defined in Chap8.tex | 3 uses, 1 definition — check ordering/spelling |
| `tab:11.1` | [Chap11.tex:76](Chap11.tex#L76) | **not defined** | table exists; `\label` missing |
| `tab:11.3` | [Chap11.tex:76](Chap11.tex#L76) | **not defined** | table exists; `\label` missing |

### `Chap6.tex:280` in detail

The source currently reads:

```latex
... transform in accordance with Eq. \eqref{4}{1}{2}.
```

This is a half-converted "Eq. 4.1(2)". `\eqref` takes one argument, so it
renders as `(??)12` — the trailing `{1}{2}` leaks into the text as a literal
"12". Should be `\eqref{eq:4:1:2}`.

The same sentence appears again at [Chap6.tex:613](Chap6.tex#L613) as the
unconverted literal text "Eq. 4.1(2)", so there is a matching pair to convert.

### Chapter 11 tables

Referenced: `tab:11.1`, `tab:11.3`, `tab:11.4`, `tab:11.6`, `tab:11.9`,
`tab:11.10`, `tab:11.13`, `tab:11.14`, `tab:11.15`.
Defined: `tab:11.2`, `tab:11.4`–`tab:11.15`.

So `tab:11.1` and `tab:11.3` are the only gaps. Separately, **9 places in
Chapter 11 write `Table~11.N` as plain text** rather than `\ref` — worth
converting for consistency, though not an error.

---

## C. Large commented-out regions

Two chapters are mostly inside `\begin{comment}` blocks and never reach the
document. Index entries were added to the live parts only (verified: zero
entries past each boundary).

| Chapter | Live | Commented out | Share | Content withheld |
|---|---|---|---|---|
| [Chap9.tex](Chap9.tex) | 1–788 | 789–1773 | **55%** | §"Numerical Values of the 6j Symbols", Tables 9.9–9.11 |
| [Chap10.tex](Chap10.tex) | 1–1374 | 1375–4129 | **67%** | Tables 10.1–10.14 |

If these come back, they need their own indexing pass. If they are staying out,
the `\section` headings inside them are dead and could be deleted rather than
commented.

---

## D. The glossary is printing as raw LaTeX source

[glossary.tex](glossary.tex) — "GLOSSARY OF SYMBOLS AND NOTATION" — has all
**271** of its content lines inside two `\begin{verbatim}` blocks. It is
`\include`d and is in the `\includeonly` list, so it is in the document, and it
currently typesets literally:

```
GLOSSARY OF SYMBOLS AND NOTATION
$\mathbf{e}_{x}, \mathbf{e}_{y}, \mathbf{e}_{z}$
$\mathbf{e}_{r}, \mathbf{e}_{\vartheta}, \mathbf{e}_{\varphi}$
```

— dollar signs, backslashes and all, rather than as typeset mathematics.
(Verified in the current PDF, not just the source.)

This looks like unconverted OCR parked in `verbatim` to keep it out of the way.
Two things to decide:

1. **Whether to convert it.** Dropping the `verbatim` wrappers would typeset the
   maths, but the file is a bare list of symbols with no descriptions — the
   right-hand column of the book's glossary appears to be missing.

2. **Whether it is still wanted at all.** The Index of Symbols built during this
   pass covers the same ground (194 entries, each with a description and page
   references) and is generated rather than hand-maintained. If the glossary
   stays, it is worth reading as a checklist of symbols the index may not yet
   reach; if it goes, `\include{glossary}` and the `\includeonly` entry come out
   with it.

---

## E. Label-naming inconsistencies

Not errors — nothing is broken by these — but they will bite anyone writing a
cross-reference by pattern.

- **Chapter 13 uses colons where the rest of the book uses dots.**
  Section labels there read `sec:13:1`, `sec:13:2.1` … rather than `sec:13.1`,
  `sec:13.2.1`. Within Chapter 13 one label is inconsistent even with its own
  siblings: [Chap13.tex:1640](Chap13.tex#L1640) is `\label{sec:13:2:8}` (colon)
  while 13.2.1–13.2.7 and 13.2.9 use `sec:13:2.N` (dot). Nothing references it.

- **Chapter 1 uses a capital.** [Chap1.tex:1](Chap1.tex#L1) is
  `\label{Chap:1}`; every other chapter is `\label{chap:N}`. Nothing references
  it, so `\ref{chap:1}` would silently fail if ever written.

- **Chapter 11 mixes reference styles for tables:** 9 places write
  `Table~11.N` as literal text rather than `\ref{tab:11.N}`.

---

## F. Indexing decisions worth a second opinion

These are judgment calls made during the indexing pass, recorded so they can be
overridden cheaply.

1. **`{j₁j₂j₃}` indexed as "triangle symbol", not "3j symbol."**
   [Chap4.tex:551](Chap4.tex#L551) calls it "the 3j-symbol", but it is the
   triad/triangle delta, not the Wigner 3j symbol that Chapter 8 defines.
   Indexing it as `3j symbol` would have merged two different objects onto one
   index line. Unified with Chapter 9's `{abc}`, so it now collects three
   chapters.

2. **Same-function/different-letter symbols were unified**, following the
   decision made for the generalized character. Affected: `χ^J_λ(ω)` (Ch2 wrote
   `χ^S_L`), `j_l(x)` (Ch4 wrote `j_λ`), `P_l^m(x)` (Ch4 wrote `P_l^λ`),
   `C_ν^λ(x)` (Ch4 wrote `C_{2J}^1`). Index keys only — the displayed maths in
   the body text is untouched everywhere.

3. **12j symbols have no symbol-index entry.** Their notation is a 3×4 array
   with an offset middle row; as a `smallmatrix` it is roughly twice the width
   of an index column. They are findable under
   `12j symbol → of the first/second kind`.

4. **`W^∥_{JM}` omitted** from the symbol index ([Chap7.tex](Chap7.tex)) while
   `W^⊥_{JM}` is present. `\|` in an index entry puts a bare `|` — makeindex's
   encapsulator character — in front of the sort key.

---

## G. Housekeeping

- **`\includeonly` in [VMK.tex](VMK.tex#L222).** Commit `8e73b04` landed during a
  build and captured a temporary value
  (`{Chap1,…,Chap10,references,glossary}`). It has been restored to
  `{Chap1,references,glossary}`, which therefore shows as an uncommitted diff.

- **Generated index files are tracked.** `.gitignore` now lists `VMK.sdx` and
  `VMK.snd` (and `VMK.idx` twice), but those files were committed in `82853e1`,
  so the ignore has no effect until `git rm --cached VMK.idx VMK.sdx VMK.snd`.

---

## Fixed during the indexing pass

For the record, these were repaired rather than logged:

| What | Where | Effect |
|---|---|---|
| `label{chap:6}` missing its backslash | [Chap6.tex:1](Chap6.tex#L1) | 7 `\ref{chap:6}` were rendering `??`; stray text in the chapter title |
| `\label{chap:12}` absent | [Chap12.tex:1](Chap12.tex#L1) | 3 `\ref{chap:12}` were rendering `??` |
| `\label{chap:9}` absent | [Chap9.tex:1](Chap9.tex#L1) | consistency; every other chapter has one |
| `\includegraphics[a;t={}]fig12_2_13b` | [Chap12.tex:437](Chap12.tex#L437) | 9 cascading LaTeX errors; Chapter 12 now builds clean |
| `\cite{ref0110}` (bib has `ref110`) | [Chap9.tex:6](Chap9.tex#L6) | 1 undefined citation |
| Two overfull index lines | Chap7, Chap8 entries | 6.4pt → 1.6pt; symbol index now clean |
