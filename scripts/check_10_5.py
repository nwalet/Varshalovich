#!/usr/bin/env python3
r"""
Checks for Section 10.5 (recursion relations for the 9j symbols) of Chapter 10,
Varshalovich, Moskalev & Khersonskii.  Each relation is tested on all valid
argument sets in a small grid where every radicand is non-negative, comparing
LHS and RHS numerically.

  eq 10.5.2/10.5.3/10.5.4   four-9j recursions (lambda=1/2)
  eq 10.5.5                 five-9j recursion (with A_q, eq 10.5.6)
  eq 10.5.7/10.5.8/10.5.9   five-9j recursions [10.5.9 tests -a fix]
  eq 10.5.10                six-9j recursion   [tests 1/3 -> 1/2 fix]
  eq 10.5.11                two equal columns
  eq 10.5.13                one argument = 0/1

Usage:  python3 check_10_5.py
"""
import math
from sympy import Rational, S, sqrt
from sympy.physics.wigner import wigner_9j

H = Rational(1, 2)


def tri(a, b, c):
    return (a >= 0 and b >= 0 and c >= 0 and abs(a - b) <= c <= a + b
            and (a + b + c) == int(a + b + c))


def valid9(v):
    a, b, c, d, e, f, g, h, j = v
    return all([tri(a, b, c), tri(d, e, f), tri(g, h, j),
                tri(a, d, g), tri(b, e, h), tri(c, f, j)])


def w9(*v):
    return wigner_9j(*v) if valid9(v) else S.Zero


def sq(*xs):
    if any(x < 0 for x in xs):
        return None
    p = S.One
    for x in xs:
        p *= x
    return sqrt(p)


def close(u, v):
    d = complex((S(u) - S(v)).evalf(30))
    return math.isfinite(d.real) and abs(d) < 1e-13


def grid(vals):
    for a in vals:
        for b in vals:
            for c in vals:
                for d in vals:
                    for e in vals:
                        for f in vals:
                            for g in vals:
                                for h in vals:
                                    for j in vals:
                                        yield (a, b, c, d, e, f, g, h, j)


VG = [1, Rational(3, 2)]


# ---- eq 10.5.2 ----
def eq2(a, b, c, d, e, f, g, h, j):
    l1 = sq(a - b + c + H, a + b + c + Rational(3, 2), c + f - j + H, c + f + j + Rational(3, 2))
    l2 = sq(-a + b + c + H, a + b - c + H, -c + f + j + H, c - f + j + H)
    r1 = sq(d + e - f + 1, -d + e + f, -a + g + d + 1, a - d + g)
    r2 = sq(d - e + f, d + e + f + 1, a + d - g, a + d + g + 1)
    if None in (l1, l2, r1, r2):
        return None
    lhs = (l1 * w9(a, b, c + H, d, e, f, g, h, j) + l2 * w9(a, b, c - H, d, e, f, g, h, j)) / (2 * c + 1)
    rhs = (r1 * w9(a - H, b, c, d + H, e, f - H, g, h, j) + r2 * w9(a - H, b, c, d - H, e, f - H, g, h, j)) / (2 * d + 1)
    return lhs, rhs


# ---- eq 10.5.3 ----
def eq3(a, b, c, d, e, f, g, h, j):
    l1 = sq(-a + b + c + H, a + b - c + H, -c + f + j + H, c - f + j + H)
    l2 = sq(a - b + c + H, a + b + c + Rational(3, 2), c + f - j + H, c + f + j + Rational(3, 2))
    r1 = sq(d - e + f + 1, d + e + f + 2, a + d - g + 1, a + d + g + 2)
    r2 = sq(-d + e + f + 1, d + e - f, -a + d + g, a - d + g + 1)
    if None in (l1, l2, r1, r2):
        return None
    lhs = (l1 * w9(a, b, c + H, d, e, f, g, h, j) + l2 * w9(a, b, c - H, d, e, f, g, h, j)) / (2 * c + 1)
    rhs = (r1 * w9(a + H, b, c, d + H, e, f + H, g, h, j) + r2 * w9(a + H, b, c, d - H, e, f + H, g, h, j)) / (2 * d + 1)
    return lhs, rhs


# ---- eq 10.5.4 ----
def eq4(a, b, c, d, e, f, g, h, j):
    l1 = sq(a - b + c + H, a + b + c + Rational(3, 2), -c + f + j + H, c - f + j + H)
    l2 = sq(-a + b + c + H, a + b - c + H, c + f - j + H, c + f + j + Rational(3, 2))
    r1 = sq(d - e + f + 1, d + e + f + 2, -a + d + g + 1, a - d + g)
    r2 = sq(-d + e + f + 1, d + e - f, a + d - g, a + d + g + 1)
    if None in (l1, l2, r1, r2):
        return None
    lhs = (l1 * w9(a, b, c + H, d, e, f, g, h, j) - l2 * w9(a, b, c - H, d, e, f, g, h, j)) / (2 * c + 1)
    rhs = (r1 * w9(a - H, b, c, d + H, e, f + H, g, h, j) - r2 * w9(a - H, b, c, d - H, e, f + H, g, h, j)) / (2 * d + 1)
    return lhs, rhs


# ---- eq 10.5.5 / 10.5.6 ----
def A(q, p, r, s, t):
    return sq(-p + r + q, p - r + q, p + r - q + 1, p + r + q + 1,
              -s + t + q, s - t + q, s + t - q + 1, s + t + q + 1)


def eq5(a, b, c, d, e, f, g, h, j):
    if c == 0 or d == 0:
        return None
    Acp, Ac = A(c + 1, a, b, f, j), A(c, a, b, f, j)
    Adp, Ad = A(d + 1, e, f, a, g), A(d, e, f, a, g)
    if None in (Acp, Ac, Adp, Ad):
        return None
    lhs = (Acp / ((c + 1) * (2 * c + 1)) * w9(a, b, c + 1, d, e, f, g, h, j)
           + Ac / (c * (2 * c + 1)) * w9(a, b, c - 1, d, e, f, g, h, j)
           - Adp / ((d + 1) * (2 * d + 1)) * w9(a, b, c, d + 1, e, f, g, h, j)
           - Ad / (d * (2 * d + 1)) * w9(a, b, c, d - 1, e, f, g, h, j))
    P = lambda x: x * (x + 1)
    rhs = ((P(a) + P(d) - P(g)) * (P(d) - P(e) + P(f)) / P(d)
           - (P(a) - P(b) + P(c)) * (P(c) + P(f) - P(j)) / P(c)) * w9(a, b, c, d, e, f, g, h, j)
    return lhs, rhs


# ---- eq 10.5.7 ----
def eq7(a, b, c, d, e, f, g, h, j):
    t1 = sq(g + h + j + 1, g + h - j, -b + e + h, b + e - h + 1)
    t2 = sq(g - h + j, -g + h + j + 1, b + e + h + 2, b - e + h + 1)
    t3 = sq(g + h + j + 2, g + h - j + 1, b + e + h + 2, b - e + h + 1)
    t4 = sq(g - h + j + 1, -g + h + j, b + e - h + 1, -b + e + h)
    D1 = sq(a + d + g + 2, a - d + g + 1)
    D2 = sq(a + d - g + 1, -a + d + g)
    rr = sq(a + d + g + 2, a - d + g + 1, a + d - g + 1, -a + d + g)
    if None in (t1, t2, t3, t4, D1, D2, rr) or rr == 0:
        return None
    lhs = (t1 / D1 * w9(a + H, b + H, c, d, e, f, g - H, h - H, j)
           + t2 / D1 * w9(a + H, b + H, c, d, e, f, g - H, h + H, j)
           + t3 / D2 * w9(a + H, b + H, c, d, e, f, g + H, h + H, j)
           - t4 / D2 * w9(a + H, b + H, c, d, e, f, g + H, h - H, j))
    rhs = (2 * g + 1) * (2 * h + 1) * sq(a + b + c + 2, a + b - c + 1) / rr * w9(a, b, c, d, e, f, g, h, j)
    return lhs, rhs


# ---- eq 10.5.8 ----
def eq8(a, b, c, d, e, f, g, h, j):
    t1 = sq(g + h + j + 1, g + h - j, b + e + h + 1, b + h - e)
    t2 = sq(-g + h + j + 1, g - h + j, -b + e + h + 1, b + e - h)
    t3 = sq(g - h + j + 1, -g + h + j, b + e + h + 1, b - e + h)
    t4 = sq(g + h + j + 2, g + h - j + 1, -b + e + h + 1, b + e - h)
    D1 = sq(a + d - g, -a + d + g + 1)
    D2 = sq(a + d + g + 1, a - d + g)
    rr = sq(a + d + g + 1, a - d + g, a + d - g, -a + d + g + 1)
    if None in (t1, t2, t3, t4, D1, D2, rr) or rr == 0:
        return None
    lhs = (t1 / D1 * w9(a - H, b - H, c, d, e, f, g - H, h - H, j)
           - t2 / D1 * w9(a - H, b - H, c, d, e, f, g - H, h + H, j)
           + t3 / D2 * w9(a - H, b - H, c, d, e, f, g + H, h - H, j)
           + t4 / D2 * w9(a - H, b - H, c, d, e, f, g + H, h + H, j))
    rhs = (2 * g + 1) * (2 * h + 1) * sq(a + b + c + 1, a + b - c) / rr * w9(a, b, c, d, e, f, g, h, j)
    return lhs, rhs


# ---- eq 10.5.9 (with -a fix in the second denominator) ----
def eq9(a, b, c, d, e, f, g, h, j):
    t1 = sq(-g + h + j + 1, g - h + j, -b + e + h + 1, b + e - h)
    t2 = sq(g + h + j + 1, g + h - j, b + e + h + 1, b - e + h)
    t3 = sq(g - h + j + 1, -g + h + j, b + e + h + 1, b - e + h)
    t4 = sq(g + h + j + 2, g + h - j + 1, -b + e + h + 1, b + e - h)
    D1 = sq(a + d + g + 2, a - d + g + 1)
    D2 = sq(a + d - g + 1, -a + d + g)
    rr = sq(a + d + g + 2, a - d + g + 1, a + d - g + 1, -a + d + g)
    if None in (t1, t2, t3, t4, D1, D2, rr) or rr == 0:
        return None
    lhs = (t1 / D1 * w9(a + H, b - H, c, d, e, f, g - H, h + H, j)
           - t2 / D1 * w9(a + H, b - H, c, d, e, f, g - H, h - H, j)
           + t3 / D2 * w9(a + H, b - H, c, d, e, f, g + H, h - H, j)
           + t4 / D2 * w9(a + H, b - H, c, d, e, f, g + H, h + H, j))
    rhs = (2 * g + 1) * (2 * h + 1) * sq(a - b + c + 1, -a + b + c) / rr * w9(a, b, c, d, e, f, g, h, j)
    return lhs, rhs


# ---- eq 10.5.10 (six-9j; tests 1/3 -> 1/2) ----
def eq10(a, b, c, d, e, f, g, h, j):
    def P2(x, y):
        return sq((x + 1) ** 2 - y ** 2) if (x + 1) ** 2 - y ** 2 >= 0 else None
    if c == 0 or f == 0 or j == 0 or (f + j - c) == 0 or (c - f + j) == 0 or (c + f - j) == 0:
        return None
    A1 = P2(a + b, c) and sq((a + b + 1) ** 2 - (c + 1) ** 2, (c + 1) ** 2 - (a - b) ** 2)
    A2 = sq((d + e + 1) ** 2 - (f + 1) ** 2, (f + 1) ** 2 - (d - e) ** 2)
    A3 = sq((g + h + 1) ** 2 - (j + 1) ** 2, (j + 1) ** 2 - (g - h) ** 2)
    B1 = sq((a + b + 1) ** 2 - c ** 2, c ** 2 - (a - b) ** 2)
    B2 = sq((d + e + 1) ** 2 - f ** 2, f ** 2 - (d - e) ** 2)
    B3 = sq((g + h + 1) ** 2 - j ** 2, j ** 2 - (g - h) ** 2)
    rf = sq(Rational(f + j - c + 1, f + j - c)) if f + j - c > 0 else None
    if None in (A1, A2, A3, B1, B2, B3, rf):
        return None
    lhs = (A1 / (2 * c + 1) * sqrt(Rational(f + j - c + 1, f + j - c)) * w9(a, b, c + 1, d, e, f, g, h, j)
           + A2 / (2 * f + 1) * sqrt(Rational(c - f + j + 1, c - f + j)) * w9(a, b, c, d, e, f + 1, g, h, j)
           + A3 / (2 * j + 1) * sqrt(Rational(c + f - j + 1, c + f - j)) * w9(a, b, c, d, e, f, g, h, j + 1))
    pre = sqrt(Rational((-c + f + j + 1) * (c - f + j + 1) * (c + f - j + 1) * (c + f + j + 2),
                        (-c + f + j) * (c - f + j) * (c + f - j) * (c + f + j + 1)))
    rhs = pre * (B1 / (2 * c + 1) * sqrt(Rational(-c + f + j, -c + f + j + 1)) * w9(a, b, c - 1, d, e, f, g, h, j)
                 + B2 / (2 * f + 1) * sqrt(Rational(c - f + j, c - f + j + 1)) * w9(a, b, c, d, e, f - 1, g, h, j)
                 + B3 / (2 * j + 1) * sqrt(Rational(c + f - j, c + f - j + 1)) * w9(a, b, c, d, e, f, g, h, j - 1))
    return lhs, rhs


# ---- eq 10.5.11 (two equal columns) ----
def eq11(a, c, d, f, g, j):
    l1 = sq(2 * a + c + 2, 2 * a - c, -c + f + j, c - f + j + 1, c + f - j + 1, c + f + j + 2)
    l2 = sq(2 * a + c + 1, 2 * a - c + 1, -c + f + j + 1, c - f + j, c + f - j, c + f + j + 1)
    r1 = sq(2 * d + f + 2, 2 * d - f, -c + f + j + 1, c - f + j, c + f - j + 1, c + f + j + 2)
    r2 = sq(2 * d + f + 1, 2 * d - f + 1, -c + f + j, c - f + j + 1, c + f - j, c + f + j + 1)
    if None in (l1, l2, r1, r2):
        return None
    lhs = (l1 * w9(a, a, c + 1, d, d, f, g, g, j) + l2 * w9(a, a, c - 1, d, d, f, g, g, j)) / (2 * c + 1)
    rhs = (r1 * w9(a, a, c, d, d, f + 1, g, g, j) + r2 * w9(a, a, c, d, d, f - 1, g, g, j)) / (2 * f + 1)
    return lhs, rhs


# ---- eq 10.5.13 (one arg = 0/1) ----
def eq13(a, b, c, d, e, g):
    if c == 0 or g == 0:
        return None
    P = lambda x: x * (x + 1)
    lhs = w9(a, b, c, d, e, c, g, g, 1)
    rhs = (P(a) + P(e) - P(d) - P(b)) / (2 * sqrt(P(c) * P(g))) * w9(a, b, c, d, e, c, g, g, 0)
    return lhs, rhs


def run():
    print("Section 10.5 recursion checks\n")
    ok = True
    for name, fn in [("eq 10.5.2", eq2), ("eq 10.5.3", eq3), ("eq 10.5.4", eq4),
                     ("eq 10.5.5", eq5), ("eq 10.5.7", eq7), ("eq 10.5.8", eq8),
                     ("eq 10.5.9 [-a fix]", eq9), ("eq 10.5.10 [1/3->1/2]", eq10)]:
        good = bad = 0
        for v in grid(VG):
            if not valid9(v):
                continue
            r = fn(*v)
            if r is None:
                continue
            good += 1
            if not close(*r):
                bad += 1
        okk = good > 0 and bad == 0
        ok &= okk
        print(f"  [{'OK  ' if okk else 'FAIL'}] {name:22s} ({good} cases, {bad} bad)")

    # eq 10.5.11 : {a a c; d d f; g g j}
    for name, fn, builder in [("eq 10.5.11", eq11, lambda a, c, d, f, g, j: (a, a, c, d, d, f, g, g, j))]:
        good = bad = 0
        for a in VG:
            for c in VG:
                for d in VG:
                    for f in VG:
                        for g in VG:
                            for j in VG:
                                if not valid9(builder(a, c, d, f, g, j)):
                                    continue
                                r = fn(a, c, d, f, g, j)
                                if r is None:
                                    continue
                                good += 1
                                if not close(*r):
                                    bad += 1
        okk = good > 0 and bad == 0
        ok &= okk
        print(f"  [{'OK  ' if okk else 'FAIL'}] {name:22s} ({good} cases, {bad} bad)")

    # eq 10.5.13 : {a b c; d e c; g g 0/1}
    good = bad = 0
    for a in VG:
        for b in VG:
            for c in VG:
                for d in VG:
                    for e in VG:
                        for g in VG:
                            if not valid9((a, b, c, d, e, c, g, g, 1)):
                                continue
                            r = eq13(a, b, c, d, e, g)
                            if r is None:
                                continue
                            good += 1
                            if not close(*r):
                                bad += 1
    okk = good > 0 and bad == 0
    ok &= okk
    print(f"  [{'OK  ' if okk else 'FAIL'}] {'eq 10.5.13':22s} ({good} cases, {bad} bad)")

    print("\n  (not checked: 10.5.1 general; 10.5.12; 10.5.14-16)")
    print("\nALL 10.5 CHECKS PASS" if ok else "\nSOME 10.5 CHECKS FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
