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


def prodf(args):
    if any(F(t) is None for t in args):
        return None
    p = S.One
    for t in args:
        p *= fac(t)
    return p


# ---- eq 10.8.10 : {d+g b c; d e f; g h c+f} ----
def eq810(b, c, d, e, f, g, h):
    a, j = d + g, c + f
    v = (a, b, c, d, e, f, g, h, j)
    if not valid9(v):
        return None
    pref = ((-1) ** (d + h - b - f) * sqrt(fac(2 * c) * fac(2 * f) * fac(2 * d) * fac(2 * g)
            / (fac(2 * c + 2 * f + 1) * fac(2 * d + 2 * g + 1))) * fac(d + e - f)
            / (D(d + g, b, c) * D(c + f, g, h))
            * fac(b + e - h) * fac(d + g + b - c) * fac(c + f + h - g)
            / (D(b, e, h) * D(d, e, f) * fac(d + e + f + 1) * fac(b + e + h + 1) * fac(d + g + b + c + 1)))
    s = S.Zero
    for x in nrange(2 * e, b + e - h, d + e - f):
        t = prodf([2 * e - x, h + b - e + x, d + g + e + c - h - x])
        u = prodf([x, b + e - h - x, d + e - f - x])
        if t is None or u is None:
            continue
        s += t / u
    return close(pref * s, w9(v))


# ---- eq 10.8.11 : {a b a+b; d e f; a+d h j} ----
def eq811(a, b, d, e, f, h, j):
    v = (a, b, a + b, d, e, f, a + d, h, j)
    if not valid9(v):
        return None
    pref = (D(a + b, f, j) * D(a + d, h, j) / (D(b, e, h) * D(d, e, f))
            * sqrt(fac(2 * b) * fac(2 * d) / (fac(2 * a + 2 * b + 1) * fac(2 * a + 2 * d + 1)))
            * fac(a + b + f + j + 1) * fac(a + d + h + j + 1) * fac(h - b + e) * fac(e - d + f)
            / (fac(f + j - a - b) * fac(h + j - a - d) * fac(b + e + h + 1) * fac(d + e + f + 1)))
    s = S.Zero
    for x in nrange(a + b + f - j, a + d + h - j):
        t = prodf([a + b + e + d - j - x, j + f - a - b + x, j + h - a - d + x])
        u = prodf([x, a + b + f - j - x, a + d + h - j - x, j + e - a - b - d + x, 2 * j + 1 + x])
        if t is None or u is None:
            continue
        s += (-1) ** x * t / u
    return close(pref * s, w9(v))


# ---- eq 10.8.12 : {a b a-b; d e f; a-d h j} (sum from x=1) ----
def eq812(a, b, d, e, f, h, j):
    if a - b < 0 or a - d < 0:
        return None
    v = (a, b, a - b, d, e, f, a - d, h, j)
    if not valid9(v):
        return None
    pref = ((-1) ** (b + f - d - h) * sqrt(fac(2 * a - 2 * b) * fac(2 * a - 2 * d) * fac(2 * b) * fac(2 * d))
            / (D(a - b, f, j) * D(a - d, h, j) * D(b, e, h) * D(d, e, f))
            * fac(b + f + j - a) * fac(d + h + j - a) * fac(h - b + e) * fac(e - d + f)
            / (fac(a - b + f + j + 1) * fac(a - d + h + j + 1) * fac(b + e + h + 1) * fac(d + e + f + 1)))
    s = S.Zero
    for x in range(1, int(2 * j) + 3):
        t = prodf([x - 1, b + e + d - a - j + x - 1, a - d + h + j - x + 1, a - b + f + j - x + 1])
        u = prodf([2 * j + 1 - x, b + f - a - j + x - 1, d + h - a - j + x - 1, a + e + j - b - d - x + 1])
        if t is None or u is None:
            continue
        s += t / u
    return close(pref * s, w9(v))


# ---- eq 10.8.13 : {a b a+b; d e f; g h a+b+f} (closed form) ----
def eq813(a, b, d, e, f, g, h):
    v = (a, b, a + b, d, e, f, g, h, a + b + f)
    if not valid9(v):
        return None
    val = ((-1) ** (a + d - g) * D(a + b + f, g, h) / (D(a, d, g) * D(b, e, h) * D(d, e, f))
           * sqrt(fac(2 * a) * fac(2 * b) * fac(2 * f) / ((2 * a + 2 * b + 1) * fac(2 * a + 2 * b + 2 * f + 1)))
           * fac(a + b + g + h + f + 1) * fac(g - a + d) * fac(e - b + h) * fac(d + e - f)
           / (fac(g + h - a - b - f) * fac(a + g + d + 1) * fac(b + e + h + 1) * fac(d + e + f + 1)))
    return close(val, w9(v))


def chk(v, val):
    if val is None:
        return None
    return close(val, w9(v))


# ---- Sec 10.8.4 three degenerate triads ----
def eq814(a, b, d, e, g, h):
    v = (a, b, a + b, d, e, d + e, g, h, g + h)
    if not valid9(v):
        return None
    return chk(v, D(a + b, d + e, g + h) / (D(a, d, g) * D(b, e, h))
               * fac(a + b + e + d + g + h + 1) / (fac(a + d + g + 1) * fac(b + e + h + 1))
               * sqrt(fac(2 * a) * fac(2 * b) * fac(2 * d) * fac(2 * e) * fac(2 * g) * fac(2 * h)
                      / (fac(2 * a + 2 * b + 1) * fac(2 * d + 2 * e + 1) * fac(2 * g + 2 * h + 1))))


def eq815(b, c, d, e, g, h):
    v = (c + b, b, c, d, e, d + e, g, h, g + h)
    if not valid9(v) or (d - b - c + g) < 0:
        return None
    return chk(v, (-1) ** (b + e - h) * D(c + b, d, g) / (D(c, d + e, g + h) * D(b, e, h))
               * fac(d + e - c + g + h) / (fac(b + e + h + 1) * fac(d - b - c + g))
               * sqrt(fac(2 * b) * fac(2 * c) * fac(2 * d) * fac(2 * e) * fac(2 * g) * fac(2 * h)
                      / (fac(2 * b + 2 * c + 1) * fac(2 * d + 2 * e + 1) * fac(2 * g + 2 * h + 1))))


def eq816(a, b, e, f, g, h):
    v = (a, b, a + b, e + f, e, f, g, h, a + b + f)
    if not valid9(v):
        return None
    return chk(v, (-1) ** (a + e + f - g) * D(a + b + f, g, h) / (D(a, g, e + f) * D(b, e, h))
               * fac(a + b + g + h + f + 1) * fac(e - b + h) / (fac(g + h - a - b - f) * fac(b + e + h + 1))
               * fac(e + g + f - a) / fac(a + e + f + g + 1)
               * sqrt(fac(2 * a) * fac(2 * b) * fac(2 * e)
                      / ((2 * a + 2 * b + 1) * fac(2 * e + 2 * f + 1) * fac(2 * a + 2 * b + 2 * f + 1))))


def eq817(a, b, c, d, f, j):
    v = (a, b, c, d, d + f, f, a + d, b + d + f, j)
    if not valid9(v):
        return None
    return chk(v, D(a + d, b + d + f, j) / (D(a, b, c) * D(c, f, j))
               * fac(a + b + 2 * d + f + j + 1) * fac(a - b + c) * fac(j + c - f)
               / (fac(a + b + c + 1) * fac(c + f + j + 1) * fac(a - b - f + j))
               * sqrt(fac(2 * a) * fac(2 * b) * fac(2 * f)
                      / (fac(2 * a + 2 * d + 1) * (2 * f + 2 * d + 1) * fac(2 * b + 2 * d + 2 * f + 1))))


# ---- Sec 10.8.5 four degenerate triads ----
def eq818(a, b, d, e, f):
    v = (a, b, a + b, d, e, f, a + d, b + e, a + b + f)
    if not valid9(v):
        return None
    return chk(v, D(a + b + f, a + d, b + e) / D(d, e, f)
               * fac(2 * a + 2 * b + d + e + f + 1) / fac(d + e + f + 1)
               * sqrt(fac(2 * e) * fac(2 * d) * fac(2 * f)
                      / ((2 * a + 2 * b + 1) * fac(2 * a + 2 * d + 1) * fac(2 * e + 2 * b + 1) * fac(2 * a + 2 * b + 2 * f + 1))))


def eq820(a, b, f, g, h):
    v = (a, b, a + b, a + g, b + h, f, g, h, a + b + f)
    if not valid9(v) or (-a - b + g + h - f) < 0:
        return None
    return chk(v, (-1) ** (2 * a) * D(a + b + f, g, h) / D(a + g, b + h, f)
               * fac(a + b + g + h - f) / fac(-a - b + g + h - f)
               * sqrt(fac(2 * f) * fac(2 * g) * fac(2 * h)
                      / ((2 * a + 2 * b + 1) * fac(2 * a + 2 * b + 2 * f + 1) * fac(2 * a + 2 * g + 1) * fac(2 * b + 2 * h + 1))))


def eq821(a, b, d, f, h):
    v = (a, b, a + b, d, b + h, f, a + d, h, a + b + f)
    if not valid9(v) or (d + h - f - b) < 0:
        return None
    return chk(v, D(h, a + d, a + b + f) / D(d, f, b + h)
               * fac(2 * a + b + d + h + f + 1) * fac(d + b + h - f) / (fac(h + d + b + f + 1) * fac(d + h - f - b))
               * sqrt(fac(2 * d) * fac(2 * f) * fac(2 * h)
                      / ((2 * a + 2 * b + 1) * fac(2 * a + 2 * d + 1) * fac(2 * b + 2 * h + 1) * fac(2 * a + 2 * b + 2 * f + 1))))


def eq822(a, b, f, g, h):
    v = (a, b, a + b, a + g, a + g + f, f, g, h, a + b + f)
    if not valid9(v):
        return None
    return chk(v, (-1) ** (2 * a) * D(a + b + f, g, h) / D(b, h, a + g + f)
               * fac(a - b + f + g + h) / fac(g + h - a - b - f)
               * sqrt(fac(2 * b) * fac(2 * g)
                      / ((2 * a + 2 * g + 1) * (2 * a + 2 * b + 1) * fac(2 * a + 2 * f + 2 * g + 1) * fac(2 * a + 2 * b + 2 * f + 1))))


def eq823(a, b, d, f, h, third):
    v = (a, b, a + b, d, d + f, f, a + d, h, a + b + f)
    if not valid9(v) or (h + d - f - b) < 0:
        return None
    return chk(v, D(a + b + f, a + d, h) / D(b, d + f, f)
               * fac(2 * a + b + d + h + f + 1) * fac(h - b + d + f) / (fac(b + d + f + h + 1) * fac(h + d - f - b))
               * sqrt(fac(2 * b) * fac(2 * d) * fac(third)
                      / ((2 * a + 2 * b + 1) * fac(2 * a + 2 * d + 1) * fac(2 * d + 2 * f + 1) * fac(2 * a + 2 * b + 2 * f + 1))))


def eq824(a, b, d, f, g):
    v = (a, b, a + b, d, d + f, f, g, b + d + f, a + b + f)
    if not valid9(v):
        return None
    return chk(v, (-1) ** (a + d - g) * D(a + b + f, g, d + f + b) / D(a, d, g)
               * fac(2 * b + 2 * f + a + d + g + 1) / fac(a + d + g + 1)
               * sqrt(fac(2 * a) * fac(2 * d)
                      / ((2 * a + 2 * b + 1) * (2 * d + 2 * f + 1) * fac(2 * a + 2 * b + 2 * f + 1) * fac(2 * d + 2 * b + 2 * f + 1))))


def eq825(a, b, d, f, h):
    v = (a, b, a + b, d, b + h, f, a + b + f + h, h, a + b + f)
    if not valid9(v) or (b + d + h - f) < 0:
        return None
    return chk(v, (-1) ** (d - b - f - h) / (D(d, f, b + h) * D(a, d, a + b + f + h))
               * fac(b + d + h - f) / (fac(2 * a + d + b + f + h + 1) * (d + b + h + f + 1))
               * sqrt(fac(2 * a) * fac(2 * f) * fac(2 * a + 2 * b + 2 * h + 2 * f + 1)
                      / ((2 * a + 2 * b + 1) * (2 * a + 2 * b + 2 * f + 1) * fac(2 * b + 2 * h + 1))))


def eq826(a, b, e, f, j, third):
    v = (a, b, a + b, e + f, e, f, a + e + f, b + e, j)
    if not valid9(v):
        return None
    return chk(v, (-1) ** (f + a + b - j) * D(j, b + e, a + e + f) / D(f, a + b, j)
               * fac(f + j + b + a + 2 * e + 1) / fac(f + a + b + j + 1)
               * sqrt(fac(2 * f) * fac(2 * b) * fac(third)
                      / ((2 * f + 2 * e + 1) * fac(2 * e + 2 * b + 1) * fac(2 * a + 2 * b + 1) * fac(2 * f + 2 * e + 2 * a + 1))))


# ---- Sec 10.8.6 five / six degenerate ----
def eq827(a, b, d, e):
    v = (a, b, a + b, d, e, d + e, a + d, b + e, a + b + d + e)
    if not valid9(v):
        return None
    return chk(v, 1 / sqrt((2 * a + 2 * b + 1) * (2 * d + 2 * e + 1) * (2 * a + 2 * d + 1) * (2 * b + 2 * e + 1)))


def eq828(a, b):
    v = (a, b, a + b, a, b, a + b, 2 * a, 2 * b, 2 * a + 2 * b)
    if not valid9(v):
        return None
    return chk(v, 1 / ((2 * a + 2 * b + 1) * sqrt((4 * a + 1) * (4 * b + 1))))


def eq829(a):
    v = (a, a, 2 * a, a, a, 2 * a, 2 * a, 2 * a, 4 * a)
    if not valid9(v):
        return None
    return chk(v, S(1) / (4 * a + 1) ** 2)


def eq830(a, b, f, g):
    v = (a, b, a + b, a + g, a + g + f, f, g, a + b + f + g, a + b + f)
    if not valid9(v):
        return None
    return chk(v, (-1) ** (2 * a) / sqrt((2 * a + 2 * g + 1) * (2 * a + 2 * b + 1) * (2 * a + 2 * b + 2 * f + 1) * (2 * a + 2 * g + 2 * f + 1)))


def eq831(a, b, f):
    v = (a, b, a + b, a + b, a + b + f, f, b, a + 2 * b + f, a + b + f)
    if not valid9(v):
        return None
    return chk(v, (-1) ** (2 * a) / ((2 * a + 2 * b + 1) * (2 * a + 2 * b + 2 * f + 1)))


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

    # eq 10.8.10 - 10.8.13 : iterate their own free-parameter grids
    P = [H, 1, Rational(3, 2), 2]

    def scan(name, fn, nvars):
        good = bad = 0
        import itertools
        for combo in itertools.product(P, repeat=nvars):
            r = fn(*combo)
            if r is None:
                continue
            good += 1
            if not r:
                bad += 1
        okk = good > 0 and bad == 0
        print(f"  [{'OK  ' if okk else 'FAIL'}] {name:24s}    ({good - bad}/{good} cases)")
        return okk

    ok &= scan("eq 10.8.10 algebraic sum", eq810, 7)
    ok &= scan("eq 10.8.11 algebraic sum", eq811, 7)
    ok &= scan("eq 10.8.12 algebraic sum", eq812, 7)
    ok &= scan("eq 10.8.13 closed form", eq813, 7)

    for name, fn, n in [("eq 10.8.14", eq814, 6), ("eq 10.8.15", eq815, 6), ("eq 10.8.16", eq816, 6),
                        ("eq 10.8.17", eq817, 6), ("eq 10.8.18", eq818, 5), ("eq 10.8.20", eq820, 5),
                        ("eq 10.8.21", eq821, 5), ("eq 10.8.22", eq822, 5), ("eq 10.8.24", eq824, 5),
                        ("eq 10.8.25", eq825, 5), ("eq 10.8.27", eq827, 4), ("eq 10.8.28", eq828, 2),
                        ("eq 10.8.29", eq829, 1), ("eq 10.8.30", eq830, 4), ("eq 10.8.31", eq831, 3)]:
        ok &= scan(name + " closed form", fn, n)

    # FLAGGED (OCR errors that could not be reconstructed numerically):
    #  10.8.19 : LHS printed as a 4x2 array, not a valid 3x3 9j.
    #  10.8.23 : an argument-dependent error remains beyond the doubled (2d)!(2d)!.
    #  10.8.26 : an argument-dependent error remains beyond the doubled (2b)!(2b)!.
    print("  [FLAG] eq 10.8.19  LHS is an OCR-mangled 4x2 array, not a valid 9j")
    print("  [FLAG] eq 10.8.23  argument-dependent OCR error (doubled (2d)! and more)")
    print("  [FLAG] eq 10.8.26  argument-dependent OCR error (doubled (2b)! and more)")

    print("\nALL CHECKED 10.8 FORMS PASS" if ok else "\nSOME 10.8 CHECKS FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
