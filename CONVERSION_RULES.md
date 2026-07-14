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

## Notes

- Rules 3-5 all use `\middle` instead of the OCR source's mismatched
  `\left`/plain-delimiter pairing, so bracket sizing scales symmetrically
  around all arguments.
- Applied so far: [Chap0.tex](Chap0.tex).
