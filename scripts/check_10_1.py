#!/usr/bin/env python3
r"""
Checks for Section 10.1 (Definition of the 9j symbols) of Chapter 10,
Varshalovich, Moskalev & Khersonskii.

Numeric (sympy):
  eq 10.1.8   9j as a sum of six Clebsch-Gordan coefficients   (definition)
  eq 10.1.9   9j orthogonality (sum over g,h)
  eq 10.1.10  9j orthogonality (sum over c,f)

Symbolic (r-symbol algebra, a..j free):
  eq 10.1.12  r_{ik}, r'_{ik} in terms of a..j
  eq 10.1.13  inverse relations (2a = r12+r13 = r'21+r'31, ...)
  eq 10.1.14  row/column sums of the r-array
  eq 10.1.15  R = a+b+c+d+e+f+g+h+j

Usage:  python3 check_10_1.py
"""
from sympy import Rational, S, sqrt, symbols, simplify
from sympy.physics.wigner import clebsch_gordan as CG, wigner_9j

H = Rational(1, 2)


def tri(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) == int(a + b + c)


def T(a, b, c):
    return 1 if tri(a, b, c) else 0


def proj(j):
    return [-j + i for i in range(int(2 * j) + 1)]


def C(a, b, c, al, be, ga):
    if a < 0 or b < 0 or c < 0 or abs(al) > a or abs(be) > b or abs(ga) > c:
        return S.Zero
    if not tri(a, b, c) or al + be != ga:
        return S.Zero
    return CG(a, b, c, al, be, ga)


def w9(a, b, c, d, e, f, g, h, j):
    if not all([tri(a, b, c), tri(d, e, f), tri(g, h, j),
                tri(a, d, g), tri(b, e, h), tri(c, f, j)]):
        return S.Zero
    return wigner_9j(a, b, c, d, e, f, g, h, j)


# --------------------------------------------------------------------------
# eq 10.1.8 : 9j = sum of six CG  (with the corrected 1/2 exponent)
# --------------------------------------------------------------------------
def eq8(j1, j2, j3, j4, j12, j34, j13, j24, j, m):
    tot = S.Zero
    for m1 in proj(j1):
        for m2 in proj(j2):
            for m3 in proj(j3):
                m4 = m - m1 - m2 - m3
                if abs(m4) > j4:
                    continue
                m12, m34, m13, m24 = m1 + m2, m3 + m4, m1 + m3, m2 + m4
                tot += (C(j1, j2, j12, m1, m2, m12) * C(j3, j4, j34, m3, m4, m34)
                        * C(j12, j34, j, m12, m34, m) * C(j1, j3, j13, m1, m3, m13)
                        * C(j2, j4, j24, m2, m4, m24) * C(j13, j24, j, m13, m24, m))
    rhs = (sqrt((2 * j12 + 1) * (2 * j13 + 1) * (2 * j24 + 1) * (2 * j34 + 1))
           * w9(j1, j2, j12, j3, j4, j34, j13, j24, j))
    return simplify(tot - rhs) == 0


# --------------------------------------------------------------------------
# eq 10.1.9 / 10.1.10 : orthogonality
# --------------------------------------------------------------------------
def rng(hi):
    return [Rational(i, 2) for i in range(0, int(2 * hi) + 1)]


def eq9(a, b, c, cp, d, e, f, fp, jj, hi=5):
    tot = S.Zero
    for g in rng(hi):
        for h in rng(hi):
            tot += ((2 * g + 1) * (2 * h + 1)
                    * w9(a, b, c, d, e, f, g, h, jj) * w9(a, b, cp, d, e, fp, g, h, jj))
    rhs = ((1 if c == cp else 0) * (1 if f == fp else 0)
           * T(a, b, c) * T(d, e, f) * T(c, f, jj) / ((2 * c + 1) * (2 * f + 1)))
    return simplify(tot - rhs) == 0


def eq10(a, b, d, e, g, gp, h, hp, jj, hi=5):
    tot = S.Zero
    for c in rng(hi):
        for f in rng(hi):
            tot += ((2 * c + 1) * (2 * f + 1)
                    * w9(a, b, c, d, e, f, g, h, jj) * w9(a, b, c, d, e, f, gp, hp, jj))
    rhs = ((1 if g == gp else 0) * (1 if h == hp else 0)
           * T(a, d, g) * T(b, e, h) * T(g, h, jj) / ((2 * g + 1) * (2 * h + 1)))
    return simplify(tot - rhs) == 0


# --------------------------------------------------------------------------
# r-symbol algebra (eq 10.1.12 - 10.1.15)
# --------------------------------------------------------------------------
def r_symbol_algebra():
    a, b, c, d, e, f, g, h, jj = symbols('a b c d e f g h jj')
    r = {
        (1, 1): -a + b + c, (1, 2): a - b + c, (1, 3): a + b - c,
        (2, 1): -d + e + f, (2, 2): d - e + f, (2, 3): d + e - f,
        (3, 1): -g + h + jj, (3, 2): g - h + jj, (3, 3): g + h - jj,
    }
    rp = {
        (1, 1): -a + d + g, (1, 2): -b + e + h, (1, 3): -c + f + jj,
        (2, 1): a - d + g, (2, 2): b - e + h, (2, 3): c - f + jj,
        (3, 1): a + d - g, (3, 2): b + e - h, (3, 3): c + f - jj,
    }
    R = a + b + c + d + e + f + g + h + jj
    out = []

    # eq 10.1.13 inverse relations
    inv = [
        (2 * a, r[1, 2] + r[1, 3]), (2 * a, rp[2, 1] + rp[3, 1]),
        (2 * d, r[2, 2] + r[2, 3]), (2 * d, rp[1, 1] + rp[3, 1]),
        (2 * g, r[3, 2] + r[3, 3]), (2 * g, rp[1, 1] + rp[2, 1]),
        (2 * b, r[1, 1] + r[1, 3]), (2 * b, rp[2, 2] + rp[3, 2]),
        (2 * e, r[2, 1] + r[2, 3]), (2 * e, rp[1, 2] + rp[3, 2]),
        (2 * h, r[3, 1] + r[3, 3]), (2 * h, rp[1, 2] + rp[2, 2]),
        (2 * c, r[1, 1] + r[1, 2]), (2 * c, rp[2, 3] + rp[3, 3]),
        (2 * f, r[2, 1] + r[2, 2]), (2 * f, rp[1, 3] + rp[3, 3]),
        (2 * jj, r[3, 1] + r[3, 2]), (2 * jj, rp[1, 3] + rp[2, 3]),
    ]
    out.append(("eq 10.1.13 inverse relations",
                all(simplify(l - rr) == 0 for l, rr in inv)))

    # eq 10.1.14 row/column sums
    sums = [
        (r[1, 1] + r[1, 2] + r[1, 3], a + b + c),
        (r[2, 1] + r[2, 2] + r[2, 3], d + e + f),
        (r[3, 1] + r[3, 2] + r[3, 3], g + h + jj),
        (r[1, 1] + r[2, 1] + r[3, 1], R - 2 * (a + d + g)),
        (r[1, 2] + r[2, 2] + r[3, 2], R - 2 * (b + e + h)),
        (r[1, 3] + r[2, 3] + r[3, 3], R - 2 * (c + f + jj)),
        (rp[1, 1] + rp[1, 2] + rp[1, 3], R - 2 * (a + b + c)),
        (rp[2, 1] + rp[2, 2] + rp[2, 3], R - 2 * (d + e + f)),
        (rp[3, 1] + rp[3, 2] + rp[3, 3], R - 2 * (g + h + jj)),
        (rp[1, 1] + rp[2, 1] + rp[3, 1], a + d + g),
        (rp[1, 2] + rp[2, 2] + rp[3, 2], b + e + h),
        (rp[1, 3] + rp[2, 3] + rp[3, 3], c + f + jj),
        (sum(r.values()), R),
        (sum(rp.values()), R),
    ]
    out.append(("eq 10.1.14 row/column sums",
                all(simplify(l - rr) == 0 for l, rr in sums)))
    return out


def run():
    print("Section 10.1 checks\n")
    ok = True

    # eq 10.1.8 : (j1,j2,j3,j4,j12,j34,j13,j24,j,m), all six triads valid
    cases8 = [(H, H, H, H, 1, 1, 1, 1, 1, 0),
              (1, 1, 1, 1, 1, 1, 1, 1, 1, 0),
              (1, H, H, 1, H, 1, Rational(3, 2), H, 1, 0)]
    r8 = all(eq8(*x) for x in cases8)
    print(f"  [{'OK  ' if r8 else 'FAIL'}] eq 10.1.8  six-CG sum = 9j  [confirms 1/2 exponent]  ({len(cases8)} cases)")
    ok &= r8

    r9 = (eq9(1, 1, 1, 1, 1, 1, 1, 1, 1) and eq9(1, 1, 1, 2, 1, 1, 1, 1, 1)
          and eq9(H, H, 1, 0, H, H, 1, 0, 1))
    print(f"  [{'OK  ' if r9 else 'FAIL'}] eq 10.1.9  orthogonality (sum g,h)")
    ok &= r9

    r10 = (eq10(1, 1, 1, 1, 1, 1, 1, 1, 1) and eq10(1, 1, 1, 1, 1, 2, 1, 1, 1)
           and eq10(H, H, H, H, 1, 0, 1, 0, 1))
    print(f"  [{'OK  ' if r10 else 'FAIL'}] eq 10.1.10 orthogonality (sum c,f)")
    ok &= r10

    print()
    for label, res in r_symbol_algebra():
        print(f"  [{'OK  ' if res else 'FAIL'}] {label}")
        ok &= res

    print("\nALL 10.1 CHECKS PASS" if ok else "\nSOME 10.1 CHECKS FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
