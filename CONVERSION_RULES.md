# Conversion rules

Running list of textual/structural conversions applied when cleaning up OCR'd
chapter files (`ChapN.tex`). Apply these to every new chapter file as it's
added. Custom commands referenced below are defined in the preamble of
[VMK.tex](VMK.tex).

## 1. Display math delimiters

Replace `$$ ... $$` with `\[ ... \]`.

```latex
% before
$$
\Psi_{\lambda\mu} = \dots
$$

% after
\[
\Psi_{\lambda\mu} = \dots
\]
```

## 2. Clebsch-Gordan coefficients → `\clebsch`

Replace `C_{a b c d}^{e f}` with `\clebsch{a}{b}{c}{d}{e}{f}`. Argument order
matches reading order: the four subscripts first, then the two superscripts.

```latex
% before
C_{j_{1} m_{1} j_{2} m_{2}}^{j m}

% after
\clebsch{j_{1}}{m_{1}}{j_{2}}{m_{2}}{j}{m}
```

Definition:

```latex
\newcommand{\clebsch}[6]{C_{#1 #2 #3 #4}^{#5 #6}}
```

## 3. Plain bra-ket → `\braket{bra}{ket}`

Replace `\left\langle ... \right| ... \rangle` (a bra and ket with a single
bar and no operator in between) with `\braket{bra}{ket}`.

```latex
% before
\left\langle a \right| b \rangle

% after
\braket{a}{b}
```

Definition:

```latex
\newcommand{\braket}[2]{\left\langle #1 \middle| #2 \right\rangle}
```

## 4. Bra-operator-ket → `\braOket{bra}{operator}{ket}`

Replace `\left\langle ... \left| ... \right| ... \rangle` (bra, an operator
sandwiched between single bars, ket) with `\braOket{bra}{operator}{ket}`.

```latex
% before
\left\langle\varepsilon^{\prime} \pi^{\prime} \alpha^{\prime} j^{\prime} m^{\prime}\right| \widehat{\mathfrak{M}}_{\lambda \mu}|\varepsilon \pi \alpha j m\rangle

% after
\braOket{\varepsilon^{\prime} \pi^{\prime} \alpha^{\prime} j^{\prime} m^{\prime}}{\widehat{\mathfrak{M}}_{\lambda \mu}}{\varepsilon \pi \alpha j m}
```

Definition:

```latex
\newcommand{\braOket}[3]{\left\langle #1 \middle| #2 \middle| #3 \right\rangle}
```

## 5. Reduced matrix element → `\braOketred{bra}{operator}{ket}`

Replace `\left\langle ... \left\| ... \right\| ... \rangle` (bra, operator
sandwiched between *double* bars, ket — the reduced-matrix-element notation)
with `\braOketred{bra}{operator}{ket}`.

```latex
% before
\left\langle\varepsilon^{\prime} \pi^{\prime} \alpha^{\prime} j^{\prime}\left\|\widehat{\mathfrak{M}}_{\lambda}\right\| \varepsilon \pi \alpha j\right\rangle

% after
\braOketred{\varepsilon^{\prime} \pi^{\prime} \alpha^{\prime} j^{\prime}}{\widehat{\mathfrak{M}}_{\lambda}}{\varepsilon \pi \alpha j}
```

Definition:

```latex
\newcommand{\braOketred}[3]{\left\langle #1 \middle\| #2 \middle\| #3 \right\rangle}
```

## 6. Heading numbering

The OCR source encodes the chapter's heading hierarchy as `\subsection*`/
`\section*` with the number baked into the title text (e.g.
`\subsection*{1.1. COORDINATE SYSTEMS. BASIS VECTORS}`). Replace these with
proper, automatically-numbered LaTeX sectioning commands, inferring depth
from how many number components the original title had:

- `\subsection*{N.M. Title}` (two-part number, one level below the chapter)
  → `\section{Title}`
- `\subsection*{N.M.K. Title}` (three-part number) → `\subsection{Title}`
- `\section*{Title}` (deepest level in the source — no numeric prefix at
  all, e.g. lettered sub-parts like `(a) ...`, or a bare title) →
  `\subsubsection{Title}`

Strip the numeric prefix (e.g. `1.1.`, `1.1.1.`, plus the space after it)
from the title text — `\section`/`\subsection` generate it automatically.
Leave lettered prefixes
like `(a)` in `\subsubsection` titles as-is, since `book` class doesn't
number `\subsubsection` by default (matching the source, which left these
deepest headings unnumbered too).

```latex
% before
\subsection*{1.1. COORDINATE SYSTEMS. BASIS VECTORS}
\subsection*{1.1.1. Cartesian Coordinate System}
\section*{(a) Description of rotations in terms of Euler angles}

% after
\section{COORDINATE SYSTEMS. BASIS VECTORS}
\subsection{Cartesian Coordinate System}
\subsubsection{(a) Description of rotations in terms of Euler angles}
```

## 7. Equation numbering

The source manually numbers equations with `\tag{n}` inside starred
environments (`equation*`, `align*`, `gather*`) or plain `\[...\]`, and
**resets `n` back to 1 at every top-level section** (i.e. at each `1.1`,
`1.2`, ... boundary — now a `\section` per rule 6). Replace this with
LaTeX's automatic numbering, reset once per *chapter* instead (so numbering
runs continuously across all of a chapter's sections), with a stable label
on every numbered equation for cross-referencing via `\eqref`.

Preamble (`VMK.tex`): `\numberwithin{equation}{chapter}` resets the
`equation` counter at each `\chapter` and renders numbers as
`chapter.equation` (e.g. `(1.5)`).

### Mechanical pass: `tools/convert_equation_numbering.py`

Run `python3 tools/convert_equation_numbering.py ChapN.tex` — it rewrites the
file in place and does the bulk of the work automatically:

- Finds every display-math block (`\[...\]`, `equation*`, `align*`,
  `gather*`) containing at least one `\tag{n}`, drops the `*` (a tagged
  `\[...\]` becomes `\begin{equation}...\end{equation}`), and replaces each
  `\tag{n}` with `\label{chap<C>:eq:<N>}`, where `<C>` is the chapter number
  (parsed from the filename) and `<N>` is a fresh sequential integer in
  document order through the whole chapter — **not** the old, per-section
  `n`, since the source resets `n` back to 1 at every top-level section and
  reusing it verbatim would produce colliding labels.
- In multi-line `align`/`gather` blocks, adds `\notag` to any row that had
  no `\tag` in the source (non-starred `align`/`gather` number every row by
  default, unlike the starred originals). Row-splitting is nesting-aware —
  it tracks `\begin{...}/\end{...}` depth so a nested `array`/`matrix`'s own
  `\\` row breaks are never mistaken for a top-level `align`/`gather` row
  break.
- Blocks with **no** `\tag` at all are left untouched (starred/unnumbered)
  — they were intentionally unnumbered in the source.
- Also writes `ChapN.eqmap.json` (gitignored, safe to delete once the manual
  pass below is done): `"<scope>,<old local n>": <new global N>` for every
  relabeled equation, where `<scope>` counts resets (0 for the first
  top-level section, 1 for the second, ...). Needed for the manual pass
  because the old numbers repeat across scopes.

### Manual pass: prose cross-references

Rewrite every prose cross-reference to a renumbered equation (`Eq. (9)`,
`Eqs. (26)-(27)`, `Equations (43) and (44)`, etc.) to use
`\eqref{chap<C>:eq:<N>}` for each number, resolving `<N>` by looking up
`"<scope>,<old n>"` in the `.eqmap.json` — determine `<scope>` from which
top-level section the reference's surrounding prose falls in (usually, but
not always, the same section the referenced equation itself was numbered
in). Watch for:

- **Explicit cross-section overrides**: e.g. `Eqs. 1.1(29)-1.1(54)` names
  section 1.1 outright — resolve against *that* section's scope, not the
  scope the sentence physically sits in.
- **Cross-chapter references**: e.g. `Eq. 2.5(4)` (chapter 2, section 5,
  equation 4) — leave as plain text; there's nothing to `\eqref` until that
  chapter exists.
- **False positives**: text that merely looks like an equation reference,
  e.g. `$a(1)$` and `$a(2)$` denoting *two different rotations*, not
  equations (1) and (2).
- **Genuine gaps in the source**: occasionally a prose reference names an
  old number that was never actually `\tag`ged anywhere in its scope (the
  source's own OCR/transcription dropped it — e.g. an unnumbered
  restatement block that should have carried its own number). Don't
  fabricate a label for these; leave the reference as plain text and flag
  it instead of guessing.

```latex
% before
\begin{equation*}
\mathbf{r}=x \mathbf{e}_{x}+y \mathbf{e}_{y}+z \mathbf{e}_{z} \tag{1}
\end{equation*}
...
as shown in Eq. (1)

% after (chapter 1, this is the chapter's first numbered equation)
\begin{equation}
\mathbf{r}=x \mathbf{e}_{x}+y \mathbf{e}_{y}+z \mathbf{e}_{z} \label{chap1:eq:1}
\end{equation}
...
as shown in \eqref{chap1:eq:1}
```

## 8. Figure numbering

The source manually numbers figures too: each `\caption{Fig. 1.N. Title}`
is preceded by `\captionsetup{labelformat=empty}` (to suppress LaTeX's own
"Figure N:" label, since the number is baked into the caption text by
hand). Unlike equations, the source's figure numbers run straight through
the chapter with no per-section resets, so this rule is simpler — no
scope-mapping needed.

Preamble (`VMK.tex`):

```latex
\numberwithin{figure}{chapter}
\renewcommand{\figurename}{Fig.}
\captionsetup{labelsep=period}
```

`\numberwithin` resets the counter per chapter (matching equations, rule
7); `\figurename` + `labelsep=period` make LaTeX's automatic label render
as `Fig. 1.1.` instead of the default `Figure 1.1:`, matching the source.

For every figure:

- Delete the `\captionsetup{labelformat=empty}` line.
- Strip the `Fig. 1.N.` (or `Fig. 1.N` — the source is inconsistent about
  the trailing period) prefix from the caption text.
- Add `\label{chap<C>:fig:<N>}` right after the `\caption{...}` (same
  `chap<C>:` namespace as equations, `:fig:` instead of `:eq:`).

Then rewrite every prose reference (`(Fig. 1.3)`, `shown in Fig. 1.5.`,
etc.) to `Fig.~\ref{chap<C>:fig:<N>}`, keeping surrounding punctuation.
This chapter's figure numbers had no gaps or cross-chapter references to
worry about (unlike rule 7's equations) — check for both anyway before
assuming it's this simple on a future chapter.

```latex
% before
\begin{figure}[h]
\begin{center}
  \includegraphics[...]{...}
\captionsetup{labelformat=empty}
\caption{Fig. 1.1. Cartesian coordinate system.}
\end{center}
\end{figure}
...
the distances between the point and coordinate planes (Fig. 1.1).

% after
\begin{figure}[h]
\begin{center}
  \includegraphics[...]{...}
\caption{Cartesian coordinate system.}
\label{chap1:fig:1}
\end{center}
\end{figure}
...
the distances between the point and coordinate planes (Fig.~\ref{chap1:fig:1}).
```

## 9. Table numbering

Same idea as rule 8, applied to `table` environments: each
`\caption{Table 1.N\\ Title}` (note the source puts the number and title on
separate lines via `\\`, unlike figures) is preceded by
`\captionsetup{labelformat=empty}`.

Preamble (`VMK.tex`): `\numberwithin{table}{chapter}` — `\tablename`
already defaults to "Table" and `\captionsetup{labelsep=period}` is shared
with rule 8 (`\captionsetup` isn't per-environment unless you scope it), so
no extra setup is needed beyond the one line.

For every table:

- Delete the `\captionsetup{labelformat=empty}` line.
- Strip the `Table 1.N\\` prefix line from the caption (the title text
  becomes the whole caption body, same as figures).
- Add `\label{chap<C>:tab:<N>}` right after the `\caption{...}`.
- Rewrite prose references (`summarized in Table 1.3.`) to
  `Table~\ref{chap<C>:tab:<N>}`.

**Number `<N>` by auto-numbering order, not the source's printed number.**
If the source has a gap — e.g. a "Table 1.2" that was never actually a
`table` environment (garbled OCR that dropped the real table and left only
a orphaned `Table 1.2\\ Title` text fragment with no `\begin{table}` around
it) — LaTeX's `\numberwithin` counter has no way to skip it, so the next
real table renders as `1.2`, not `1.3`. Name that table's label
`chap<C>:tab:2` (matching what it actually renders as), not `chap<C>:tab:3`
(matching the source's now-stale printed number) — same principle as rule
7's equations, where labels track the new auto-numbered sequence, never the
old manual one. Leave the orphaned fragment itself as plain, unconverted
text (nothing to number) and flag it same as rule 7's "genuine gaps" case.

```latex
% before
\begin{table}[h]
\begin{center}
\captionsetup{labelformat=empty}
\caption{Table 1.1\\
Matrix form of the transformations for vector components in different bases.}
...
\end{table}
...
Table 1.2\\
Matrices of transformations between cartesian, spherical contravariant and polar components of vectors.
...
\begin{table}[h]
\begin{center}
\captionsetup{labelformat=empty}
\caption{Table 1.3\\
Differential operations.}
...
\end{table}
...
The above equations are summarized in Table 1.3.\\

% after
\begin{table}[h]
\begin{center}
\caption{Matrix form of the transformations for vector components in different bases.}
\label{chap1:tab:1}
...
\end{table}
...
Table 1.2\\
Matrices of transformations between cartesian, spherical contravariant and polar components of vectors.
...
\begin{table}[h]
\begin{center}
\caption{Differential operations.}
\label{chap1:tab:2}
...
\end{table}
...
The above equations are summarized in Table~\ref{chap1:tab:2}.\\
```

## 10. Vectors, matrices, and named operators

The source is inconsistent about the semantic weight of `\mathbf` and
`\boldsymbol`, and uses `\operatorname` for some names amsmath already
provides. Replace with two new commands (`VMK.tex`):

```latex
\newcommand{\vect}[1]{\mathbf{#1}}
\newcommand{\mat}[1]{\boldsymbol{#1}}
```

- `\mathbf{X}` → `\vect{X}` — this covers essentially every occurrence
  (single-letter vectors: `e`, `A`, `B`, `r`, `n`, ...). If a single
  `\mathbf{...}` ever wraps *two* symbols at once (an OCR artifact, e.g.
  `\mathbf{n r}` meaning the dot product `\mathbf{n}\mathbf{r}`), split it
  into two separate `\vect{...}` calls first — `\vect` should always take
  exactly one vector's name, not a whole sub-expression.
- `\boldsymbol{X}` → `\mat{X}`, but **only when `X` is actually a matrix**.
  Don't convert blindly: check each occurrence's meaning first, because the
  source also uses `\boldsymbol` for other bold quantities that aren't
  matrices — e.g. `\boldsymbol{\nabla}` (the nabla vector-operator, typeset
  bold because the source bolds vector-valued quantities throughout — this
  should become `\vect{\nabla}`, not `\mat{\nabla}`) or
  `\boldsymbol{\varepsilon}_{ikl}` (the rank-3 Levi-Civita tensor, which
  fits neither `\vect` nor `\mat` — leave it as plain `\boldsymbol`). Also
  watch for OCR noise like `\boldsymbol{\prime}` (a bold prime mark, which
  should just become a plain `\prime}`/`'` — bolding a prime is never
  meaningful).
- `\operatorname{name}` → the built-in command, **only if amsmath already
  defines `name` as a named operator with identical rendering**. In this
  chapter that's just `\operatorname{det}` → `\det`. Don't convert names
  amsmath doesn't define (`curl`, `div`, `grad` have no amsmath equivalent)
  or names where amsmath's version renders differently (`\Re`/`\Im` are
  Fraktur symbols ℜ/ℑ, not the upright-text "Re"/"Im" that
  `\operatorname{Re}`/`\operatorname{Im}` produce — swapping would be a
  visual regression, not a preserving rewrite).

```latex
% before
\mathbf{r}=x \mathbf{e}_{x}+y \mathbf{e}_{y}+z \mathbf{e}_{z}
\boldsymbol{X}^{\boldsymbol{\prime}}
\boldsymbol{\nabla} \times \mathbf{A}
-\operatorname{det} X=\mathbf{r}^{2}

% after
\vect{r}=x \vect{e}_{x}+y \vect{e}_{y}+z \vect{e}_{z}
\mat{X}^{\prime}
\vect{\nabla} \times \vect{A}
-\det X=\vect{r}^{2}
```

## 11. No blank lines around display math

The OCR source surrounds nearly every display-math block (`equation`,
`align`, `gather`, `\[...\]`) with a blank line both before and after
(often two). A blank line in LaTeX source starts a new paragraph
(`\par`), which isn't wanted here — in most cases the equation is part of
the sentence's flow, not a paragraph boundary, and the extra `\par` risks
an unwanted first-line indent or spacing artifact on the text that
follows. Delete the blank line(s) immediately before the block's
`\begin{...}`/`\[` and immediately after its `\end{...}`/`\]`, so the
equation sits flush against the surrounding prose in the source. Leave
blank lines that separate actual paragraphs (i.e., that aren't adjacent to
an equation) untouched.

```latex
% before
of a point $\vect{r}$ may be written as


\begin{equation}
\vect{r}=x \vect{e}_{x}+y \vect{e}_{y}+z \vect{e}_{z} \label{chap1:eq:1}
\end{equation}


\begin{figure}[tbh]

% after
of a point $\vect{r}$ may be written as
\begin{equation}
\vect{r}=x \vect{e}_{x}+y \vect{e}_{y}+z \vect{e}_{z} \label{chap1:eq:1}
\end{equation}
\begin{figure}[tbh]
```

## Notes

- **Compilation now requires `xelatex`** (or another Unicode engine), not
  `pdflatex` — `VMK.tex`'s preamble uses `fontspec`/`polyglossia` for font
  handling. Compile with `xelatex VMK.tex` (run 2-3 times to resolve all
  cross-references), not `pdflatex`.
- Rules 3-5 all use `\middle` instead of the OCR source's mismatched
  `\left`/plain-delimiter pairing, so bracket sizing scales symmetrically
  around all arguments.
- Applied so far: [Chap0.tex](Chap0.tex) (rules 1-5), [Chap1.tex](Chap1.tex),
  [Chap2.tex](Chap2.tex), [Chap3.tex](Chap3.tex), [Chap4.tex](Chap4.tex) and
  [Chap5.tex](Chap5.tex) (all rules — Chap2.tex and Chap3.tex have no
  figures/tables, so rules 8-9
  didn't apply there).
- Chap5.tex had two more variants of the "heading collapsed to plain text"
  bug (see Chap4.tex's notes above): a numbered `N.N.N Title` heading with
  no `\subsection*` wrapper at all (same fix — wrap it), and one case where
  the wrapper was missing *and* the title text ran straight into the
  next paragraph's prose with no separator (`5.8.1 Action of ... on
  $Y_{lm}$ Spherical components of the operator...` — split into a heading
  line and a new paragraph by eye, there's no mechanical marker for where
  the title ends and the prose begins).
- Also found `\[...\]` wrapped around a `\begin{align*}`/`\begin{gather*}`
  (nesting a numbered amsmath environment inside plain display math) —
  amsmath rejects this outright (`Erroneous nesting of equation
  structures`). Source of the bug: `$$`/`\[...\]` was OCR'd as a wrapper
  around what should have been the *only* delimiter; strip the outer
  `\[`/`\]` and keep just the inner environment. `grep -n '\\\[\s*$'`
  followed by checking the next non-blank line for `\begin{align*|gather*|equation*}`
  catches this before it surfaces as a compile error.
- Another `\tag{$\{n\}$}` (same garbled-tag bug as Chap4) — this time
  inside a `\[...\]` (single-equation) block, not an `align`/`gather` row,
  so the whole block silently stayed unconverted rather than corrupting a
  row. Checking for this pattern (`grep -n '\\tag{\$'`) *before* running
  `convert_equation_numbering.py`, not after, meant no post-hoc renumbering
  shift was needed this time.
- Rule 3's plain `\braket` appeared for the first time (previous chapters
  only had rule 4's `\braOket`) — same source pattern, `\langle a \mid
  b\rangle`, just using `\mid` instead of a bare `|` before the source's
  `\rangle`; both mean the same thing here.
- **Chap4.tex was the most structurally damaged chapter so far** — beyond
  the routine rule application, it needed:
  - **~35 more missed `(a)/(b)/(c)` subsubsection headings**, flattened to
    plain text like the rest, but *not* caught by the mechanical rule 6
    regex because they never had a numeric `N.M` prefix to match on in the
    first place (they're the deepest level, same as rule 6's `\section*`
    case, just never wrapped in *any* sectioning command). Distinguishing
    these from genuine prose that happens to start "(a) ...", "(b) ..." is
    a judgment call: a real heading is short (no full descriptive
    sentence), stands in its own paragraph, and is followed immediately by
    dedicated equation content; a prose enumeration is a multi-sentence
    paragraph that continues discursively. When in doubt, check how the
    *other* letters in the same `(a)(b)(c)...` sequence read — treat
    siblings consistently rather than one becoming a heading and another
    staying prose.
  - **Two headings trapped *inside* a `gather`/`align` body** as a
    `\text{...}` or bare line marked `\notag` (e.g. `\text { (d) Periodicity
    } \\` sitting between two tagged rows) — the OCR lost the environment
    boundary entirely, not just the sectioning command. Fix by splitting
    the math environment in two at that point and lifting the heading out.
  - **Three tables missing their entire `\begin{table}`/`\caption` wrapper**
    (just a bare `\begin{center}...\end{center}` with a plain-text title
    line like `Table 4.20.\\` in front of it) and **one table with a
    caption that belonged to a different table entirely** (copy-paste/OCR
    mixup — the caption text didn't match the tabulated data at all;
    resolved by using the orphaned title text sitting in front of the
    table, which *did* match). Six orphaned "what this table shows"
    formula lines (sitting as bare text before a table with no title of
    its own) needed merging into that table's caption.
  - **`(Cont.)` continuation tables** (a table's data split across two
    `\begin{table}` blocks because it didn't fit on one page) need to
    *not* get their own `\label`/number — `\caption` always steps the
    counter even under `labelformat=empty`, so giving each continuation a
    normal rule-9 treatment silently shifts every later table's number.
    Instead, drop `\caption` entirely for the continuation and write
    `\textbf{Table~\ref{<base's label>} (continued)}` by hand.
  - **Two equations whose `\tag{n}` was itself garbled** —
    `\tag{$\{53\}$}` instead of `\tag{53}` (stray OCR-inserted math-mode
    and escaped braces around the number). `tools/convert_equation_numbering.py`'s
    tag regex only matches pure digits, so these were silently invisible
    to it: one sat inside an already-tagged `gather` (which the script
    still processed, just skipping this one row as if untagged — no
    error, wrong-looking output) and the other was the *only* tag in its
    `equation*` block, so the script found zero tags and left the whole
    block starred and unconverted (also no error). Both only surfaced by
    grepping `\\tag\{` for leftovers *after* the script ran and expecting
    zero matches. Fixing them after the main numbering pass already ran
    means hand-computing the shift: every label at or after the insertion
    point moves up by one (by two, after the second fix) — safe to do
    because the prose `\eqref` pass hadn't run yet, so there were no
    stale cross-references to chase down too. Do the `\tag\{` leftover
    check *before* the prose pass for this reason.
  - `tools/convert_equation_numbering.py` itself had a real bug, now
    fixed: its row-splitter tracked `\begin{...}/\end{...}` nesting but
    not generic brace groups, so `\sum_{\substack{a \\ b}}` inside a
    `gather`/`align` row had its internal `\\` (meant to line-wrap the
    subscript) mistaken for a top-level row break, corrupting the
    surrounding row and — in one case — causing amsmath's own "Multiple
    \label's" error downstream. The splitter now tracks two independent
    depths: named environments, and raw (non-escaped) `{`/`}` brace
    nesting. Escaped `\{`/`\}` (the literal brace glyphs used in e.g. 6j-
    symbol notation, `\left\{j_1 j_2 j_3\right\}`) are explicitly excluded
    from the brace count. Run `grep -n '\\substack{' ChapN.tex` before
    trusting a chapter's equation numbering if it wasn't already checked.
- Chap3.tex needed two more one-off source fixes before the general rules
  applied cleanly: a `\subsection*{...}` heading whose title text itself
  contained a literal `\\` line break (merge onto one line before running
  rule 6's regex — it's anchored on end-of-line); and rule 4's bra-ket
  conversion caught a genuinely mismatched-delimiter OCR bug (`\left(` ...
  `\right\rangle` — parenthesis mixed with angle bracket), which was fixed
  in the same pass since applying `\braOket` there required first figuring
  out which side of the equation the mismatch was on.
- Rule 2's Clebsch-Gordan conversion needs whitespace between adjacent
  index digits — Chap3.tex had `C_{J M_{1} J M_{2}}^{00}` (superscript
  `j=0,m=0` written as one glued token `00`) and `C_{10 l_{i-1} 0}^{...}`
  (same issue, `10` meaning `j=1,m=0`). The brace-aware converter splits on
  whitespace and expects exactly 4 subscript / 2 superscript tokens, so it
  loudly flags these (`UNEXPECTED TOKEN COUNT`) rather than mis-converting
  — insert the missing space and rerun.
- VMK.tex's own preamble/numbering choices (e.g. `\numberwithin{equation}`
  scoped to `section` rather than `chapter`) have since been adjusted
  directly by the project owner; treat the live file as authoritative over
  this document's rule 7 example numbers if they ever disagree.
- Three known gaps left as plain text in Chap1.tex, per rule 7's / rule 9's
  "genuine gaps"/orphaned-fragment cases: "A detailed form of (25) is..."
  and "...analogous to (25)-(28)." (section 1.1 and 1.2 respectively — the
  source never tagged the equations these numbers would refer to), and the
  orphaned "Table 1.2" text fragment (never a real `table` environment, so
  nothing to number — see rule 9).
- `\boldsymbol{\varepsilon}_{ikl}` (Levi-Civita tensor) and
  `\boldsymbol{\Phi}` (a scalar field, oddly bolded once in a table row
  label) in Chap1.tex were deliberately left as plain `\boldsymbol` per
  rule 10 — neither fits `\vect` (rank 1) or `\mat` (2D matrix). Same in
  Chap2.tex for `\boldsymbol{\Theta}` (a bolded scalar angle).
- Rule 10 also applies to bolded numerals, not just letters — Chap2.tex had
  `\widehat{\mathbf{L}}^{\mathbf{2}}` (a bold "2" exponent, OCR noise like
  rule 10's `\boldsymbol{\prime}` case): drop the bold entirely rather than
  wrapping a digit in `\vect` or `\mat`.
- Environment-boundary detection (`tools/convert_equation_numbering.py`)
  requires `\begin{...}`/`\end{...}` to be alone on their line. Chap2.tex
  had one `\end{align*} \quad(i, k, l=x, y, z)` with trailing content stuck
  on the same line (OCR noise) — the script silently merged it with the
  next matching `\end{align*}` much later in the file, corrupting both
  blocks' numbering, with no error at conversion time (it only surfaced as
  a `\begin{align} on input line N ended by \end{align*}` compile error).
  Grep for `\\end\{(equation\*|align\*|gather\*)\}\s*\S` (and the
  `\begin` equivalent with leading content) on every new chapter *before*
  running the script, and move any trailing content onto the equation's
  last content line first.
