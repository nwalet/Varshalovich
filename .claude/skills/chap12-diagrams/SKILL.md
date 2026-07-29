---
name: chap12-diagrams
description: Convert the scanned recoupling/coupling diagrams in the VMK book's Chapter 12 (Sec 12.1 3jm-sums, Sec 12.2 6j/9j-sums) into line-only TikZ that replaces the OCR image crops in Chap12.tex. Use when the task is redrawing, inserting, or fixing any Chap12 diagram (6j-diamonds, 9j-hexagons, 12j octagons, thetas, loops, delta-rects), or when the user references t12dev.tex, diagrams.sty macros, image crops (p7..p42 or 39a53bec-* hashes), or "the next diagram".
---

# Chapter 12 scanned-diagram → TikZ conversion

Redraw the book's hand-drawn recoupling diagrams as clean, line-only TikZ.
Every diagram is an equation `LHS = RHS` between nj-symbol graphs. Reproduce
the **book's picture** (crop is authoritative for layout), but derive ambiguous
interior labels from the **equation** via the triad rule below.

## Golden rules
- **Lines + labels only.** No arrows, no node dots, no +/- signs.
- **Line weight encodes summation.** A momentum SUMMED in the equation
  (printed as a CAPITAL letter: X, Y, Z, X_1...) is an **ultra-thick** line
  (`dlu`). Everything else is **thick** (`dl`).
- Match the crop's shape, orientation, and label placement; preserve symmetry.
- Compile clean under `xelatex VMK.tex` after every insert.

## Workflow (dev loop)
1. **Inspect** the crop: `Read images/<name>.(jpg|png)`. Identify each piece
   (see Vocabulary) and read the equation in Chap12.tex for exact labels.
2. **Build** in `t12dev.tex` (a `standalone` scratch file). One `\[...\]` per
   diagram; use `\begin{scope}[shift={(x,y)}]` (or `xshift=`) per piece.
3. **Render**: `xelatex -interaction=nonstopmode t12dev.tex` then
   `pdftoppm -png -r 150 t12dev.pdf /tmp/out`; Read the PNG and compare to crop.
4. **Review**: the user co-edits t12dev in their IDE (especially interior
   labels). Re-read t12dev before inserting — they often tweak it.
5. **Insert**: replace the matching `\includegraphics{...}` (or its
   `\begin{center}...\end{center}` block) in Chap12.tex with the `\[...\]`.
   Map crop→equation by section/equation order.
6. **Verify**: `xelatex VMK.tex`, confirm `grep -c "^!"` is 0, and that the
   remaining `grep -c includegraphics Chap12.tex` count dropped by one.
7. **Commit** at each subsection boundary (message: which eqs, which shapes).
   Remote is `Varshalovich`; push with `git push Varsalovich HEAD`.

Working preference (strategy 2): build fast, insert clean-looking matches
directly, flag only genuinely uncertain interiors for the user.

## Macros (all in diagrams.sty; do NOT redefine locally)
Styles: `dl` (thick j-line), `dlu` (ultra-thick, summed), `wb` (white-bg label).
- `\dloop{x}{y}{lab}` — closed circle = Pi_lab factor.
- `\cthetav[midstyle]{x}{y}{left}{mid}{right}` — theta as circle + VERTICAL diameter.
- `\ctheta[midstyle]{x}{y}{top}{mid}{bot}` — theta as circle + HORIZONTAL diameter.
- `\drect{x}{y}{top}{bot}` / `\drects{x}{y}{top}{bot}{left}{right}` — delta rectangle.
- `\dsixjsq[N]{x}{y}{tl}{tr}{bl}{br}{hd}{vd}` — 6j diamond/square, 2 diagonals.
  `N` = numeric bold selector: 0 both thin, 1 vertical(vd) bold, 2 horizontal(hd)
  bold, 3 both bold. Exposes `qt,qr,qb,ql` (top/right/bottom/left corners) so you
  can overdraw a single EDGE bold, e.g. `\draw[dlu](ql)--(qb);` bolds the BL edge.
- `\dsq{x}{y}{top}{bot}{left}{right}` — axis-aligned square; exposes
  `sqTL,sqTR,sqBL,sqBR`. Draw diagonals yourself:
  `\draw[dl](sqTL)--(sqBR)node[wb,pos=0.3]{$..$};` (node-on-draw = label near line).
- `\dninehex{x}{y}{eTL}{eL}{eBL}{eBR}{eR}{eTR}` — POINTY-top 9j hexagon (6 edges
  clockwise from upper-left); exposes `hT,hUL,hLL,hB,hLR,hUR`. Draw the 3 centre
  diameters yourself (bold the summed one): `\draw[dlu](hT)--(hB)node[wb,pos=0.27]{$Z$};`
- `\dhexflat{x}{y}{top}{ur}{lr}{bot}{ll}{ul}` — FLAT-top 9j hexagon; exposes
  `v0..v5` (v0=right, CCW). Diameters: v0-v3, v1-v4, v2-v5.
- `\dsixjtri`, `\dtheta`, `\dlens`, `\dthetav`, `\dlensv` — older Sec 12.1 shapes.
- 6j-TRIANGLE (apex-up/left/right) is usually inlined: `\coordinate(c)at(0,0);`
  + 3 vertices + 3 outer edges + 3 spokes `\draw[dl](c)--(V)node[wb,pos=0.5]{..}`.

## Vocabulary (which symbol → which shape)
- Pi_a factor → **loop** (circle, `\dloop`).
- delta / triangle {a b c} → **theta** (circle+diameter, `\cthetav`/`\ctheta`)
  or a small **rect** (`\drect`) when it's a Kronecker delta between two lines.
- 6j `{. . .;. . .}` → **diamond** (`\dsixjsq`) or **square-with-diagonals**
  (`\dsq` + 2 diagonals) or **triangle-with-centre** (inline).
- 9j → **hexagon** (pointy `\dninehex` OR flat `\dhexflat`) with 3 centre diameters.
- 12j **first kind** (`\twelvejI`) → **octagon-WHEEL**: 8 outer edges + 4 full
  diameters through ONE centre point.
- 12j **second kind** (`\twelvejII`) → **octagon-LATTICE**: 8 outer edges + a
  central `#` (2 verticals + 2 horizontals) forming a small central square.
  (Octagon coords pattern: Ta/Tb top, Ra/Rb right, Qa/Qb bottom, La/Lb left, at
  (+/-0.6,+/-1.4)/(+/-1.4,+/-0.6); wheel uses 8 vertices at 45-deg steps radius 1.6.)

## Triad rule — resolve ambiguous INTERIOR labels from the equation
Each nj vertex is a TRIAD of three j-lines that MEET at a point.
- **9j** `{a b c / d e f / g h i}`: every ROW and every COLUMN is a triad
  ({a b c},{d e f},{g h i},{a d g},{b e h},{c f i}). In the hexagon the 6 outer
  vertices realise these (2 adjacent edges + 1 diameter each). K_{3,3} form: the
  outer 6-cycle = the (i,j) entries R1-C1-R2-C2-R3-C3; the 3 diameters = the
  remaining R1-C2, R2-C3, R3-C1 entries.
- **6j** `{a b c / d e f}`: four triads = {a b c}, {a e f}, {d b f}, {d e c}.
See memory `chap12-nj-vertex-triads` and `chap12-diagram-conventions`.

## Bold a single outer edge
Macros draw all outer edges thin. To bold ONE summed edge, overdraw between the
exposed corner/vertex coords AFTER the macro, e.g. for a diamond BL edge
`\draw[dlu](ql)--(qb);`; for a flat-hex top edge `\draw[dlu](v2)--(v1);`.

## Status / handoff
DONE — Chapter 12 is fully converted: `grep -c includegraphics Chap12.tex` == 0,
`xelatex VMK.tex` clean. Sec 12.1 (eq 12.1.2–29) and Sec 12.2 (eq 30–63) are all
line-only TikZ. If reopening: the diagrams may still get cosmetic label/geometry
tweaks on review — same dev loop applies. This skill + the two memory files are
the reference for any future nj-diagram work.
