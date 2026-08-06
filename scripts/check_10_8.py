#!/usr/bin/env python3
r"""
Checks for Section 10.8 (explicit 9j forms at special relations between
arguments -- degenerate triads) of Chapter 10, Varshalovich, Moskalev &
Khersonskii.  Closed forms are compared with sympy's wigner_9j.

This file is built up incrementally.  Covered so far:
  eq 10.8.2   one degenerate triad (double algebraic sum)
  eq 10.8.3   the same, as a single CG sum

Usage:  python3 check_10_8.py
"""
import math
from sympy import Rational, S, sqrt, factorial as fac
from sympy.physics.wigner import clebsch_gordan as CG, wigner_9j

H = Rational(1, 2)


def tri(a, b, c):
    return (a >= 0 and b >= 0 and c >= 0 and abs(a - b) <= c <= a + b
            and (a + b + c) == int(a + b + c))


def valid9(v):
    a, b, c, d, e, f, g, h, j = v
    return all([tri(a, b, c), tri(d, e, f), tri(g, h, j),
                tri(a, d, g), tri(b, e, h), tri(c, f, j)])


def w9(v):
    return wigner_9j(*v) if valid9(v) else S.Zero


def D(a, b, c):
    return sqrt(fac(a + b - c) * fac(a - b + c) * fac(-a + b + c) / fac(a + b + c + 1))


def C(a, b, c, al, be, ga):
    if abs(al) > a or abs(be) > b or abs(ga) > c or al + be != ga or not tri(a, b, c):
        return S.Zero
    return CG(a, b, c, al, be, ga)


def F(x):
    return fac(x) if (x >= 0 and x == int(x)) else None


def close(u, v):
    d = complex((S(u) - S(v)).evalf(30))
    return math.isfinite(d.real) and abs(d) < 1e-12


def nrange(*mx):
    return range(0, int(max(mx)) + 2)


# ---- eq 10.8.3 : one degenerate triad, single CG sum ----
def eq83(a, b, c, d, e, f, h, j):
    g = a + d
    if not valid9((a, b, c, d, e, f, g, h, j)):
        return None
    pref = (sqrt(fac(a + d + j - h) * fac(h + j + a + d + 1)) / (D(a, b, c) * D(d, e, f) * D(c, f, j))
            * sqrt(fac(2 * a) * fac(2 * d) / (fac(2 * a + 2 * d + 1) * (2 * h + 1)))
            * fac(c + b - a) * fac(f + e - d) * fac(f + c - j)
            / (fac(a + b + c + 1) * fac(f + e + d + 1) * fac(f + c + j + 1)))
    s = S.Zero
    for x in nrange(2 * c, f + c - j):
        d1, d2 = f + c - j - x, x
        if d1 < 0 or F(f + j - c + x) is None or F(2 * c - x) is None:
            continue
        base = (-1) ** x * fac(f + j - c + x) * fac(2 * c - x) / (fac(x) * fac(f + c - j - x))
        rn, rd = b + a - c + x, e + d + c - j - x
        dn, dd = b + c - a - x, e + j - c - d + x
        if min(rn, rd, dn, dd) < 0:
            continue
        rad = sqrt(fac(rn) * fac(rd) / (fac(dn) * fac(dd)))
        cg = C(b, e, h, c - a - x, j - d - c + x, j - a - d)
        s += base * rad * cg
    return close(pref * s, w9((a, b, c, d, e, f, g, h, j)))


# ---- eq 10.8.2 : one degenerate triad, double algebraic sum (corrected) ----
# OCR fixes: prefactor (b+e+h+1)!(d+e+f+1)! (were b+c+h, d+c+f) and the sum
# denominator (c-e-a+h-x+y) (was the corrupt "c-c-a").
def eq82(a, b, c, d, e, f, h, j):
    g = a + d
    if not valid9((a, b, c, d, e, f, g, h, j)):
        return None
    pref = (D(h, j, a + d) / (D(a, b, c) * D(b, e, h) * D(d, e, f) * D(c, f, j))
            * sqrt(fac(2 * a) * fac(2 * d) / fac(2 * a + 2 * d + 1))
            * fac(b + c - a) * fac(b + e - h) * fac(e + f - d) * fac(c + f - j) * fac(a + d + h + j + 1)
            / (fac(a + b + c + 1) * fac(b + e + h + 1) * fac(d + e + f + 1) * fac(c + f + j + 1)))
    s = S.Zero
    for x in nrange(2 * c):
        for y in nrange(2 * e):
            args_num = [2 * c - x, 2 * e - y, j + f - c + x, h + b - e + y]
            args_den = [x, y, c + f - j - x, b + e - h - y, c - e - a + h - x + y, e - c - d + j + x - y]
            if any(F(t) is None for t in args_num + args_den):
                continue
            num = S.One
            for t in args_num:
                num *= fac(t)
            den = S.One
            for t in args_den:
                den *= fac(t)
            s += (-1) ** (x + y) * num / den
    return close(pref * s, w9((a, b, c, d, e, f, g, h, j)))


CASES = [(1, 1, 1, 1, 1, 1, 1, 1), (1, 1, 2, 1, 1, 1, 1, 2), (1, 2, 2, 1, 1, 1, 2, 1),
         (H, 1, H, 1, 1, 1, Rational(3, 2), 1), (1, 1, 1, 1, 2, 2, 1, 2),
         (Rational(3, 2), 1, H, 1, 1, 1, 1, Rational(3, 2)), (1, 2, 1, 1, 1, 2, 2, 1)]


# ---- eq 10.8.8 : two degenerate triads {a b c; d e f; a+d b+e j}, algebraic sum ----
def eq88(a, b, c, d, e, f, j):
    g, h = a + d, b + e
    if not valid9((a, b, c, d, e, f, g, h, j)):
        return None
    pref = (D(a + d, b + e, j) / (D(a, b, c) * D(d, e, f) * D(c, f, j))
            * sqrt(fac(2 * a) * fac(2 * b) * fac(2 * d) * fac(2 * e)
                   / (fac(2 * a + 2 * b + 1) * fac(2 * d + 2 * e + 1)))
            * fac(a + b + d + e + j + 1) * fac(a - b + c) * fac(d - e + f) * fac(c + f - j)
            / (fac(a + b + c + 1) * fac(d + e + f + 1) * fac(c + f + j + 1)))
    s = S.Zero
    for z in nrange(2 * f, c + f - j):
        ad = [c + f - j - z, d + f - e - z, a - b - f + j + z]
        if F(2 * f - z) is None or F(j + c - f + z) is None or any(F(t) is None for t in ad):
            continue
        s += ((-1) ** z * fac(2 * f - z) * fac(j + c - f + z)
              / (fac(z) * fac(c + f - j - z) * fac(d + f - e - z) * fac(a - b - f + j + z)))
    return close(pref * s, w9((a, b, c, d, e, f, g, h, j)))


# ---- eq 10.8.9 : same, single CG (confirms the /->f fix) ----
def eq89(a, b, c, d, e, f, j):
    g, h = a + d, b + e
    if not valid9((a, b, c, d, e, f, g, h, j)):
        return None
    rn = [2 * a, 2 * b, 2 * d, 2 * e, a + b + d + e + j + 1, a + d + e + b - j]
    rd = [2 * a + 2 * d + 1, 2 * b + 2 * e + 1, a + b + c + 1, a + b - c, d + e + f + 1, d + e - f]
    if any(F(t) is None for t in rn + rd):
        return None
    num = S.One
    for t in rn:
        num *= fac(t)
    den = (2 * j + 1)
    for t in rd:
        den *= fac(t)
    val = sqrt(num / den) * C(c, f, j, a - b, d - e, a - b + d - e)
    return close(val, w9((a, b, c, d, e, f, g, h, j)))


def run():
    print("Section 10.8 checks (incremental)\n")
    ok = True

    r83 = [eq83(*x) for x in CASES]
    r83 = [r for r in r83 if r is not None]
    o83 = all(r83) and len(r83) > 0
    print(f"  [{'OK  ' if o83 else 'FAIL'}] eq 10.8.3  CG-sum form          ({sum(r83)}/{len(r83)} cases)")
    ok &= o83

    rr = [r for r in (eq82(*x) for x in CASES) if r is not None]
    o82 = all(rr) and len(rr) > 0
    ok &= o82
    print(f"  [{'OK  ' if o82 else 'FAIL'}] eq 10.8.2  double-sum [corrected] ({sum(rr)}/{len(rr)} cases)")

    for name, fn, args in [("eq 10.8.8  algebraic sum", eq88, [(a, b, c, d, e, f, j) for (a, b, c, d, e, f, h, j) in CASES]),
                           ("eq 10.8.9  CG [/->f fix]", eq89, [(a, b, c, d, e, f, j) for (a, b, c, d, e, f, h, j) in CASES])]:
        rs = [r for r in (fn(*x) for x in args) if r is not None]
        oo = all(rs) and len(rs) > 0
        ok &= oo
        print(f"  [{'OK  ' if oo else 'FAIL'}] {name:24s}    ({sum(rs)}/{len(rs)} cases)")

    print("\n(remaining 10.8.10-10.8.31 pending)")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
