#!/usr/bin/env python3
r"""
Checks for Section 10.2 of Chapter 10 (explicit forms of the 9j symbols),
Varshalovich, Moskalev & Khersonskii.

  eq 10.2.15  9j as a sum of six CG   [tests the epsilon-corrected form]
  eq 10.2.16  9j as a sum of six CG   [tests the epsilon-corrected form]
  eq 10.2.17  9j as a sum of six 3jm
  eq 10.2.20  9j as a sum of three 6j
  eq 10.2.21  9j as a sum of three Racah W

Not checked (structural only): 10.2.1 (quadruple algebraic sum, heavily OCR-
damaged), 10.2.2-10.2.14 (Wu 15-variable sums), 10.2.18/10.2.19.

Usage:  python3 check_10_2.py
"""
import math
from sympy import Rational, S, sqrt, factorial as fac
from sympy.physics.wigner import clebsch_gordan as CG, wigner_3j, wigner_6j, wigner_9j, racah

H = Rational(1, 2)


def tri(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) == int(a + b + c)


def proj(j):
    return [-j + i for i in range(int(2 * j) + 1)]


def C(a, b, c, al, be, ga):
    if a < 0 or b < 0 or c < 0 or abs(al) > a or abs(be) > b or abs(ga) > c:
        return S.Zero
    if not tri(a, b, c) or al + be != ga:
        return S.Zero
    return CG(a, b, c, al, be, ga)


def w3(a, b, c, al, be, ga):
    if abs(al) > a or abs(be) > b or abs(ga) > c or al + be + ga != 0 or not tri(a, b, c):
        return S.Zero
    return wigner_3j(a, b, c, al, be, ga)


def w6(a, b, c, d, e, f):
    if not (tri(a, b, c) and tri(a, e, f) and tri(d, b, f) and tri(d, e, c)):
        return S.Zero
    return wigner_6j(a, b, c, d, e, f)


def w9(a, b, c, d, e, f, g, h, j):
    if not all([tri(a, b, c), tri(d, e, f), tri(g, h, j),
                tri(a, d, g), tri(b, e, h), tri(c, f, j)]):
        return S.Zero
    return wigner_9j(a, b, c, d, e, f, g, h, j)


def Wv(p, q, r, s, t, u):
    if not (tri(p, q, t) and tri(s, r, t) and tri(p, r, u) and tri(s, q, u)):
        return S.Zero
    return racah(p, q, r, s, t, u)


def valid9(a, b, c, d, e, f, g, h, j):
    return all([tri(a, b, c), tri(d, e, f), tri(g, h, j),
                tri(a, d, g), tri(b, e, h), tri(c, f, j)])


def close(u, v):
    d = complex((S(u) - S(v)).evalf(30))
    return math.isfinite(d.real) and abs(d) < 1e-13


# ---- eq 10.2.15 : six CG (epsilon-corrected) ----
def eq15(a, b, c, d, e, f, g, h, j):
    tot = S.Zero
    for al in proj(a):
        for be in proj(b):
            for de in proj(d):
                for ep in proj(e):
                    ga, ph, eta, mu = al + be, de + ep, al + de, be + ep
                    nu = ga + ph
                    tot += (C(a, b, c, al, be, ga) * C(d, e, f, de, ep, ph)
                            * C(c, f, j, ga, ph, nu) * C(a, d, g, al, de, eta)
                            * C(b, e, h, be, ep, mu) * C(g, h, j, eta, mu, nu))
    pref = sqrt(S(1) / ((2 * c + 1) * (2 * f + 1) * (2 * g + 1) * (2 * h + 1))) / (2 * j + 1)
    return close(pref * tot, w9(a, b, c, d, e, f, g, h, j))


# ---- eq 10.2.16 : six CG (epsilon-corrected) ----
def eq16(a, b, c, d, e, f, g, h, j):
    tot = S.Zero
    for ga in proj(c):
        for be in proj(b):
            for de in proj(d):
                for mu in proj(h):
                    al = ga + be
                    eta = al - de
                    ep = be + mu
                    ph = ep - de
                    nu = mu + eta
                    tot += (C(c, b, a, ga, be, al) * C(g, d, a, eta, de, al)
                            * C(b, h, e, be, mu, ep) * C(d, f, e, de, ph, ep)
                            * C(h, g, j, mu, eta, nu) * C(f, c, j, ph, ga, nu))
    pref = S.NegativeOne ** (2 * (c + g)) / ((2 * a + 1) * (2 * e + 1) * (2 * j + 1))
    return close(pref * tot, w9(a, b, c, d, e, f, g, h, j))


# ---- eq 10.2.17 : six 3jm ----
def eq17(a, b, c, d, e, f, g, h, j):
    tot = S.Zero
    for al in proj(a):
        for be in proj(b):
            for de in proj(d):
                for ep in proj(e):
                    ga, ph, eta, mu = -al - be, -de - ep, -al - de, -be - ep
                    nu = -eta - mu
                    tot += (w3(a, b, c, al, be, ga) * w3(d, e, f, de, ep, ph)
                            * w3(g, h, j, eta, mu, nu) * w3(a, d, g, al, de, eta)
                            * w3(b, e, h, be, ep, mu) * w3(c, f, j, ga, ph, nu))
    return close(tot, w9(a, b, c, d, e, f, g, h, j))


# ---- eq 10.2.20 : three 6j ----
def eq20(a, b, c, d, e, f, g, h, j):
    tot = sum((-1) ** (2 * x) * (2 * x + 1)
              * w6(a, b, c, f, j, x) * w6(d, e, f, b, x, h) * w6(g, h, j, x, a, d)
              for x in [Rational(i, 2) for i in range(0, int(2 * (a + b + c + d + e + f + g + h + j)) + 1)])
    return close(tot, w9(a, b, c, d, e, f, g, h, j))


# ---- eq 10.2.21 : three Racah W ----
def eq21(a, b, c, d, e, f, g, h, j):
    tot = sum((2 * x + 1) * Wv(a, e, c, h, x, b) * Wv(a, e, g, f, x, d) * Wv(c, h, f, g, x, j)
              for x in [Rational(i, 2) for i in range(0, int(2 * (a + b + c + d + e + f + g + h + j)) + 1)])
    return close(tot, w9(a, b, c, d, e, f, g, h, j))


def D(a, b, c):
    return sqrt(fac(a + b - c) * fac(a - b + c) * fac(-a + b + c) / fac(a + b + c + 1))


# ---- eq 10.2.1 : quadruple algebraic sum (last factorial in DENOMINATOR) ----
def eq1(a, b, c, d, e, f, g, h, j):
    from functools import reduce
    pref = (D(a, b, c) * D(d, e, f) * D(b, e, h) * D(g, h, j) / (D(a, d, g) * D(c, f, j))
            * fac(a + d - g) * fac(c + f - j) * fac(g + h + j + 1)
            / (fac(a + d + g + 1) * fac(a - b + c) * fac(-a + b + c) * fac(d - e + f)
               * fac(-d + e + f) * fac(b - e + h) * fac(-b + e + h)))
    N = int(2 * (a + b + c + d + e + f + g + h + j)) + 2
    tot = S.Zero
    for x in range(N):
        for y in range(N):
            for z in range(N):
                for t in range(N):
                    nums = [2 * a - x, 2 * b - y, 2 * d - z, 2 * e - t, -a + b + c + x,
                            -b + e + h + y, -d + e + f + z, b - e + g - j + t,
                            -a - e + f + g + x + t, c - d + e + j + z - t]
                    dens = [x, y, z, t, a + b - c - x, b + e - h - y, d + e - f - z,
                            b + e - g + j - t, a + d - g - x - z, e - b + h + y - t,
                            -d + e + f + z - t, b - e + g - j - y + t, -a + c - e + g - j + x + t,
                            -a + c - d + f + g + j + 1 + x + z]
                    if any(v < 0 for v in nums + dens):
                        continue
                    num = reduce(lambda p, v: p * fac(v), nums, S(1))
                    den = reduce(lambda p, v: p * fac(v), dens, S(1))
                    tot += (-1) ** (x + y + z + t) * num / den
    return close(pref * (-1) ** (a - c + e - g + j) * tot, w9(a, b, c, d, e, f, g, h, j))


CASES = [(H, H, 1, H, H, 1, 1, 1, 1), (1, 1, 1, 1, 1, 1, 1, 1, 1),
         (1, H, H, 1, H, H, 1, 1, 1), (Rational(3, 2), 1, H, 1, 1, 1, H, 1, Rational(3, 2)),
         (1, 1, 2, 1, 1, 1, 1, 1, 1), (H, 1, H, 1, H, H, H, H, 1)]


def run():
    print("Section 10.2 checks\n")
    ok = True
    for name, fn in [("eq 10.2.1  quad sum [denom-corrected]", eq1),
                     ("eq 10.2.15 six-CG [eps-corrected]", eq15),
                     ("eq 10.2.16 six-CG [eps-corrected]", eq16),
                     ("eq 10.2.17 six-3jm", eq17),
                     ("eq 10.2.20 three-6j", eq20),
                     ("eq 10.2.21 three-Racah", eq21)]:
        good = bad = 0
        for x in CASES:
            if not valid9(*x):
                continue
            good += 1
            if not fn(*x):
                bad += 1
        okk = good > 0 and bad == 0
        ok &= okk
        print(f"  [{'OK  ' if okk else 'FAIL'}] {name:36s} ({good} cases, {bad} bad)")
    print("\n  (not checked: 10.2.1 quad sum; 10.2.2-14 Wu; 10.2.18/19)")
    print("\nALL 10.2 CHECKS PASS" if ok else "\nSOME 10.2 CHECKS FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
