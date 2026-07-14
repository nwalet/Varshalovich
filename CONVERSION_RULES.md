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
  [Chap2.tex](Chap2.tex), [Chap3.tex](Chap3.tex), [Chap4.tex](Chap4.tex),
  [Chap5.tex](Chap5.tex), [Chap6.tex](Chap6.tex), [Chap7.tex](Chap7.tex) and
  [Chap8.tex](Chap8.tex) (all rules — Chap2.tex, Chap3.tex and Chap6.tex have
  no figures/tables, so rules 8-9 didn't apply there; Chap7.tex has no
  figures, so only rule 9 applied there).
- **Chap8.tex's back half (numerical-value tables, roughly the last third of
  the chapter) is wrapped in `\begin{comment}...\end{comment}`** by the
  project owner (predates this pass, `comment` package already loaded in
  `VMK.tex`) to keep compile times down — per explicit instruction, none of
  rules 1-11 were applied inside that block; it was left byte-for-byte
  untouched. Everything below refers only to the unwrapped front matter
  (defintions, symmetry/recursion relations, and the algebraic — not
  numerical — coefficient tables, Sec. 8.1-8.13's prose plus Tables
  8.1-8.10).
- **Chap8.tex was by far the most heavily corrupted chapter converted so
  far**, an order of magnitude more than Chap4.tex's previous high-water
  mark, because its content is almost entirely dense algebraic formulas
  (few complete sentences to anchor OCR against) using terse single-letter
  momenta (`a, b, c` / `α, β, γ`) instead of the `j_1, m_1, j_2, m_2` style
  seen elsewhere, so a lot of inter-symbol whitespace load-bears for
  tokenization in a way OCR handles poorly:
  - Rule 2's Clebsch-Gordan converter needed a genuinely iterative
    repair strategy, not just the single-space-insertion precedent: of
    ~473 `C_{...}^{...}` instances found, 313 were already well-formed, and
    a further ~90 were recovered in three escalating passes — inserting a
    space before a glued trailing `-\gamma`/`-\beta`/`-b`/etc. projection
    (the same bug as prior chapters, just far more common here); an
    unambiguous multi-split search (trying every way to split a deficient
    token list down to the required 4 subscript + 2 superscript count, only
    accepting the fix when exactly one such split existed); and, once the
    user confirmed mid-session that a bare `00` inside a Clebsch coefficient
    always means two separate zero arguments glued together (not a single
    "00"), a pass that also generalized to `060` as `0 b 0` (`b`, the second
    momentum letter, OCR'd as the digit `6`) and `+10`/`+20`-style trailing
    zeros as a glued `+1 0`/`+2 0` — each confirmed against the specific
    formula's own surrounding context (e.g. a $P_a P_b P_c$ product on the
    same line proving a `b` term must be present) rather than applied
    blindly. That left 68 genuinely unconvertible instances — mostly
    Sec. 8.5.1(i)-(k) and 8.6.3's `±`/`∓`
    template formulas (which encode two or four sibling equations at once
    through a single argument list, so no single 4+2 split is
    "the" correct one), Sec. 8.7.3-8.7.6's sums-of-three/four-CG-coefficient
    identities (corrupted well beyond spacing — stray/wrong Greek letters,
    swapped sub/superscripts), and Sec. 8.11's explicit table of *other
    authors'* differently-shaped notations (`C^{j_1j_2j_3}_{m_1m_2m}`-style,
    never meant to fit `\clebsch`'s 4+2 template at all) — left as literal
    `C_{...}^{...}$ text and not force-converted, per the existing
    "flag rather than guess" precedent (see Chap6.tex's note above).
  - Rule 9 hit its worst case yet: **two entire tables (8.5 and 8.8) had
    lost their whole `\begin{table}`/`\caption` wrapper**, not just the
    `\captionsetup` line — bare `\begin{center}\begin{tabular}` preceded by
    only a stray `Table 8.N.` text fragment (8.5) or a `\section*{Table
    8.8.}` heading plus a loose `$$formula$$` (8.8), each with its own
    correctly-wrapped `(Cont.)` half sitting right after. A **third table
    (8.10) had no textual trace at all** — not even an orphaned `Table
    8.10.` fragment — its number and formula were recoverable only because
    the formula itself survived, misplaced as a `\multicolumn` header row
    *inside* the table body. All three were reconstructed as proper
    `\begin{table}[h]\begin{center}\caption{...}\label{...}` blocks. Once
    every one of Tables 8.1-8.10 had a real `\caption`, this chapter also
    needed the *other* half of Chap4.tex's continuation precedent applied
    at scale for the first time: seven more `(Cont.)` tables (8.3, 8.4, 8.6,
    8.7, and two each for 8.9 and 8.10) each already had their own
    `\captionsetup{labelformat=empty}`+`\caption{Table 8.N. (Cont.)}` pair,
    which would have silently stepped `\numberwithin{table}`'s counter an
    extra time per continuation — all ten were rewritten to
    `\textbf{Table~\ref{chap8:tab:N} (continued)}` with no `\caption` of
    their own, so the auto-numbered sequence lands on exactly Table 8.1
    through Table 8.10 with no gaps, matching the source's own numbers.
  - Also found, and fixed the same way as rule 2's digit-glue precedent:
    **~104 instances of `]^{1/2}` (the square-root exponent that appears on
    almost every formula in this chapter) missing their leading `1`** —
    OCR'd as `]^{/2}` or `]^{/3}` — plus 11 more misread as `]^{7/2}`, 3 as
    `]^{1/3}`, 1 as `]^{1/8}`, and 1 as `]^{//2}`. All were confirmed as the
    same "square root of a factorial ratio" notation (never actually a
    7th/cube/8th root anywhere in this book) by their identical structural
    position closing a `\left[\frac{...!...!}{...!}\right]` group, and
    fixed with a global substitution rather than case-by-case, since the
    pattern and fix were unambiguous throughout.
  - Rule 6 was simpler than Chap7.tex's: **every** `8.N` and `8.N.M` heading
    was uniformly `\subsection*` regardless of its true depth (no
    inconsistent mix), so reclassifying purely by the number of
    dot-separated components in the printed title (not by which wrapper
    command the OCR happened to use) was a clean, single pass. The deepest
    unnumbered level had one genuine "sibling inconsistency" case, handled
    the same way as Chap7.tex's `(a) Coordinate inversion` /
    `(b) Rotations...` pair: **`(a) Permutations of columns` /
    `(b) Permutations of rows` / `(c) Transposition`** (Sec. 8.4.1) and
    **`(a) Rotation of the coordinate system.` / `(b) Inversion...` /
    `(c) Time reversal.`** (Sec. 8.4.5) each had only their `(b)`/`(c)`
    member surviving OCR with a `\section*{...}` wrapper; the other
    siblings, structurally identical short noun-phrase titles, were
    unwrapped plain text and were wrapped to match. Elsewhere in the
    chapter, bare `(a)/(b)/(c)/(d)` markers introducing a list of
    conditions or a group of unnumbered formula variants (dozens of them)
    were left as plain text, per Chap2.tex's precedent — this book's
    `(a)/(b)/(c)` markers are usually paragraph/equation-group labels, not
    their own headings, and none of those had a wrapped sibling as evidence
    otherwise.
  - Six more headings (`8.2.1`, `8.7.3`, `8.7.4`, `8.7.5`, `8.7.6`, `8.9.1`)
    had the familiar "collapsed to plain text" bug (no sectioning command
    at all); one heading (`8.11`) had its title text itself split across
    two source lines with the second half left as an unwrapped continuation
    line — merged onto one line before running rule 6, same fix as
    Chap3.tex's embedded-`\\`-in-a-heading case.
  - Two equations (`\tag{18a}`, `\tag{18b}` in Sec. 8.9.3, a physically
    genuine split into an $S^2\ge 0$ case and an $S^2\le 0$ case sharing one
    conceptual equation number) used a letter-suffixed tag that
    `tools/convert_equation_numbering.py`'s pure-digit regex doesn't match
    by design, so both blocks were silently left starred by the mechanical
    pass. Converted by hand into two separate `\begin{gather}...\end{gather}`
    blocks with sequential labels, which meant manually renumbering the
    two labels the automated pass had already assigned immediately after
    them (a `+2` shift, done before the prose `\eqref` pass so there was
    nothing stale to chase afterward) — same principle as Chap4.tex's
    hand-fixed `\tag{$\{53\}$}` gaps.
  - Rule 10 found the chapter's only vector-valued `\mathbf` usage
    (`\mathbf{j}`/`\mathbf{j}_1`/`\mathbf{j}_2}`, the angular-momentum
    vectors in Sec. 8.1.1) alongside a handful of **bold *scalars* used
    purely for emphasis**, not vectors or matrices: `\mathbf{J}` (and once
    the glued `\mathbf{1 3}`) labelling the perimeter value `J=a+b+c` in
    Sec. 8.10's list of vanishing `3jm` symbols (inconsistently mixed with
    plain, unbolded `J=11` a few lines later — left the unbolded instances
    alone, matching Chap7.tex's ε/ϵ precedent of only touching what's
    already marked bold); a bolded division slash,
    `Arguments Change by $\mathbf{1} \boldsymbol{/} \mathbf{2}$` as a
    heading title (dropped all three bold wrappers, rendering plain
    `$1/2$`); and `\boldsymbol{c}` bolding a table column header letter in
    four different tables (dropped, same as Chap7.tex's bold table-header
    `J` precedent). `\operatorname{Re}` (used correctly, Fraktur exception
    from rule 10 doesn't apply since this is the upright-text Re) was left
    untouched.
  - One stray pair of running header/page-number lines
    (`Clebsch-Gordan Coefficients and 3jm Symbols NNN`, OCR'd from the
    original PDF's header, three occurrences) had leaked into the body text
    between tables in the algebraic-tables section — deleted as pure OCR
    noise, not part of the actual content.
  - Two `\langle...\}`/`\{...\rangle` delimiter-mismatch bugs (an angle
    bracket OCR'd where the matching close was already a curly brace, or
    vice versa) inside otherwise-fine table-cell formulas were fixed the
    same way as Chap7.tex's `\left.\left\lvert` bug — by comparing against
    the sibling cell's matching, correctly-delimited `\{...\}` pair in the
    same table row/column.
- **Chap7.tex's heading wrapper was inconsistent with its title's own printed
  depth in two genuine source spots**: three of its four top-level `N.N`
  sections (`7.2`, `7.3`, `7.4`) were OCR'd as `\subsection*` instead of
  `\section*`, and `7.1.1. Definition` (a three-part number) was OCR'd as
  `\section*` instead of `\subsection*`. Rather than trust the source
  command, the conversion script classified every heading purely by how many
  dot-separated components its own printed number had (2 → `\section`, 3 →
  `\subsection`, none → `\subsubsection`), which handled both mismatches
  automatically. That approach then introduced its *own* bug: matching
  `\d+(?:\.\d+)*\.\s*` against a heading whose number had *no* trailing
  period before the title (e.g. `7.1.2 Components...`, vs. the more common
  `7.1.1. Definition`) let the greedy match backtrack to a shorter 2-part
  number, misclassifying 9 ordinary `\subsection*{7.N.M ...}` headings
  (`7.1.2`, `7.1.4`, `7.1.9`, `7.2.2`, `7.2.4`, `7.2.6`, `7.2.8`, `7.3.2`,
  `7.3.4` — all correctly wrapped in the source) as `\section` with a stray
  leading digit stuck in the title (`\section{2 Components of Tensor
  Spherical Harmonics}`). Caught by grepping `^\\section\{[0-9]` after the
  pass and reclassifying those 9 back to `\subsection` by hand; a correct
  regex would require the optional trailing period to be followed by
  whitespace, not swallow a period that's actually a separator before the
  next digit group.
- Chap7.tex is built from three top-level sections (`7.1`, `7.2`, `7.3`) plus
  a short unnumbered appendix-style section (`7.4`, cross-referencing other
  authors' notations) whose own three subsections (`Tensor spherical
  harmonics`, `Spinor spherical harmonics`, `Vector spherical harmonics`)
  are bare `\section*{Title}` with no numeric prefix at all — the deepest-level
  bare-title case from rule 6, converted to `\subsubsection`.
- **A whole heading (`7.3.10. Clebsch-Gordan Series`) had lost its
  `\subsection*{...}` wrapper entirely**, appearing as bare text — same
  "collapsed to plain text" bug documented for Chap4.tex/Chap5.tex. Only
  one instance this time, caught by eyeballing the full heading list rather
  than a mechanical grep (it has no lettered or numeric-but-unwrapped marker
  distinguishing it from ordinary prose at a glance).
- **The `(a) Coordinate inversion` / `(b) Rotations of coordinate system(s)`
  pair of `\subsubsection`-level headings appears three times** (once each
  in `7.1.4`, `7.2.4`, `7.3.4`, all three sections titled "Transformations of
  Coordinate Systems"). Only one of the six occurrences (`7.1.4`'s "(b)")
  survived OCR with its `\section*{...}` wrapper intact; the other five were
  flattened to plain text. Used the one surviving wrapped instance as
  evidence that this specific `(a)`/`(b)` pair is a genuine heading pair
  (short noun-phrase titles, matching Chap1.tex's precedent), then wrapped
  the other five to match. By contrast, **the many *other* `(a)`/`(b)`/`(c)`
  markers in this chapter — both bare (introducing a group of equations with
  no title text at all) and the ones followed by a full lead-in sentence
  ("(a) Spherical contravariant components ... may be written ... as")** were
  all left as plain text, per Chap2.tex's established precedent that this
  book's `(a)/(b)/(c)` markers usually label paragraphs/equation groups
  within a subsection rather than naming their own headings. No surviving
  wrapped example existed anywhere in the chapter for these, unlike the
  Coordinate-inversion pair.
- Rule 1 (`$$...$$` → `\[...\]`) needs to run independently of rule 7's
  equation-numbering script — `tools/convert_equation_numbering.py` only
  touches blocks containing a `\tag`, so Chap7.tex's seven `$$...$$` blocks
  (all untagged continuations of a numbered `align`/`gather`, alternating
  with proper environments to work around an amsmath multi-page-alignment
  limitation) were silently skipped by every earlier check and only
  surfaced by a final `grep -n '\$\$'` sweep *after* the numbering pass and
  the "no blank lines" pass had already run. Converted and had rule 11
  (blank-line stripping) re-applied to just those seven blocks by hand,
  since the earlier blanket rule-11 pass predated their conversion to
  `\[...\]` and couldn't have matched them.
- Rule 2's Clebsch-Gordan converter hit its "insert a missing space, reject
  on unexpected token count" precedent (Chap3.tex) far more than any prior
  chapter — 13 of Chap7.tex's 56 `C_{...}^{...}` instances had glued
  subscript tokens needing a space inserted (e.g. `C_{L 010}^{L'0}` →
  `C_{L 0 1 0}^{L'0}` for $L,m{=}0,S{=}1,\sigma{=}0$; `C_{10 L 0}^{J0}` →
  `C_{1 0 L 0}^{J0}$ for $j_1{=}1,m_1{=}0$; `C_{L M+\mu S-\mu}^{JM}` →
  `C_{L M+\mu S -\mu}^{JM}$ with the sign glued onto $S$). One pair
  (`C_{J M 2 n 0}^{JM}`, `C_{J 02 n 0}^{J0}`) needed the opposite fix —
  *joining* an OCR-split `2 n` into the single symbolic token `2n` (an even
  integer, not two separate quantum numbers) — which a token-count check
  alone can't catch, since `J 02 n 0` already had 4 raw whitespace-separated
  groups, just wrongly divided (`02`/`n` instead of `0`/`2n`). Verified by
  hand against the surrounding formula's physical meaning (a $J$-with-$2n$
  coupling in a Legendre-polynomial expansion) before fixing.
- Two more OCR content bugs, beyond rule 2/6's routine cases, needed
  hand reconstruction (both flagged and fixed, not guessed blind, since the
  correct form was directly evidenced by sibling equations in the same
  block): `\left.\left\lvert\, \Omega...\right.\right]^{\mp\frac12}` (a
  mismatched-delimiter bug — displays a stray `|` where the two neighboring
  rows in the same `align*` clearly show `\left[\Omega...\right]^{\mp\frac12}`
  should be there) in Sec. 7.2.2, and `\mathbf{n}\cdot[\mathbf{a} * \mathbf{a}]`
  (missing `\times` and a dropped superscript star) in two of six sibling
  equations in Sec. 7.3.9 whose other four correctly show
  `\mathbf{n}\cdot[\mathbf{a}^{*} \times \mathbf{a}]`.
- One equation's argument was garbled into an entire run-on sentence stuffed
  inside a single spurious `\mathbf{...}`:
  `$\mathbf{F ( \mathbf { r } _ { 2 } ) \Phi ( \mathbf { r } _ { 1 } , \mathbf
  { r } _ { 2 } ) \text { in a series of } \mathbf { Y } _ { J M } ^ { L }
  ( \vartheta _ { 1 } , \varphi _ { 1 } ) \text { has the form }}$` in Sec.
  7.3.14 — reconstructed as ordinary prose with two separate inline-math
  spans, `$\mathbf{F}(\mathbf{r}_2)\Phi(\mathbf{r}_1,\mathbf{r}_2)$ in a
  series of $\mathbf{Y}_{JM}^L(\vartheta_1,\varphi_1)$ has the form`, since
  every word and symbol needed was already present, just wrapped wrong.
- One prose cross-reference, `Eqs. 7.1(27)7.1(29)` in Sec. 7.3.7 (missing
  its range dash — presumably `7.1(27)-7.1(29)`), names old-numbering
  eq. 7.1(29) — but section 7.1's own `\tag` sequence genuinely skips 29
  (jumps `\tag{28}` straight to `\tag{30}` inside one `gather*`, an
  original-numbering gap, not an OCR artifact, confirmed by the surrounding
  block's own structure). Left as plain, unconverted text and flagged per
  rule 7's "genuine gaps" case, rather than guessing whether the intended
  target was `7.1(28)` (an OCR digit swap) or something else.
- Rule 10: `\boldsymbol{\nabla}` (appears ~14 times, mostly in
  `\operatorname{grad}/\operatorname{div}/\operatorname{curl}` identities of
  Sec. 7.3.6) all became `\vect{\nabla}` per the rule's own worked example.
  `\boldsymbol{\epsilon}`/`\boldsymbol{\varepsilon}` (a photon polarization
  vector in Sec. 7.3.14) became `\vect{\epsilon}`/`\vect{\varepsilon}` since
  it's unambiguously a vector — even though the *same* quantity appears
  bold in some equations and plain (unbolded) `\varepsilon(\mathbf{k})` in
  others nearby; per this chapter's convention of only touching what's
  already marked bold, the plain instances were left alone rather than
  bolding them to match. A one-off `\boldsymbol{J}` (bolding the scalar
  quantum number $J$ in a table's column header) was dropped to plain `J`,
  same as Chap2.tex's bold-numeral-exponent precedent — it fits neither
  `\vect` (not rank 1) nor `\mat` (not a matrix).
- `\operatorname{grad}`, `\operatorname{div}`, `\operatorname{curl}` (Sec.
  7.3.6) and `\operatorname{rank}` (Sec. 7.1.1) were all left as
  `\operatorname{...}` per rule 10 — amsmath has no built-in equivalent for
  any of them, matching the rule's existing `curl`/`div`/`grad` exception.
  No `\operatorname{det}` occurred in this chapter.
- Chap6.tex reused rule 7's "checks before running the script" habit to
  good effect: `grep`ing for trailing content after `\end{...}` (the
  Chap2.tex boundary-detection bug) caught 2 instances before the script
  ran, both fixed the same way (move the trailing text onto the last
  content line) with no post-hoc cleanup needed.
- Rule 4's `\braOket` conversion needed a targeted regex this time instead
  of hand-picked replacements — Chap6.tex has ~26 instances of `\langle+|
  OP|+\rangle`-style matrix elements (bra/ket are always literally `+` or
  `-`, a spin-up/spin-down label) in dense tabular/aligned blocks. Watch
  for a *different*, unrelated use of `\langle...\rangle` in the same
  chapter: `\langle X \rangle` (one pair, no `|`) is expectation-value/
  ensemble-average notation, not an inner product — it doesn't match
  `\braket` (needs a bra *and* ket) or `\braOket` (needs a middle
  operator) and should be left alone.
- One `C_{...}^{...}` was left unconverted because the fix needed to
  reconstruct actual missing content (a dropped `t_{LM}` factor and a
  superscript/subscript merged into the same wrong text), not just
  whitespace — beyond what rule 2's "insert a missing space" precedent
  covers. Flagged rather than guessed; see Chap6.tex around `\tag{50}`'s
  replacement.
- Found `\operatorname{n}` where every other occurrence of that exact unit
  vector in the same chapter is `\mathbf{n}(\vartheta,\varphi)` — an OCR
  command mix-up, not a real operator. Fixed to `\mathbf{n}` (then
  `\vect{n}` under rule 10) by cross-checking against the chapter's own
  consistent usage elsewhere, the same kind of evidence that's justified
  fixes in earlier chapters.
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
