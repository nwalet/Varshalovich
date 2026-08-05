#!/usr/bin/env python3
r"""
Numerical check of the explicit special-value formulas for the Clebsch-Gordan
coefficients listed in Section 8.5 of Varshalovich, Moskalev & Khersonskii,
"Quantum Theory of Angular Momentum".

For every formula we draw random angular momenta / projections that satisfy the
special constraint of that formula (e.g. c = a+b, or alpha = a, or
alpha = beta = gamma = 0), evaluate the book's closed-form right-hand side, and
compare it with the value returned by sympy's clebsch_gordan.  Draws that fall
outside a formula's domain (invalid projection, negative factorial argument,
...) are rejected, so only genuine instances are tested.

Convention (book == sympy):
    C_{a alpha, b beta}^{c gamma} = <a alpha, b beta | c gamma>
        == sympy.physics.wigner.clebsch_gordan(a, b, c, alpha, beta, gamma)

Covered
    8.5.1  Special values of momenta a,b,c   eq. (8.5.1)-(8.5.23)
    8.5.2  Special values of projections     eq. (8.5.32)-(8.5.48)

Not covered (defined through auxiliary functions C,D,E,F that still carry OCR
artifacts -- flag for manual repair, as was done for the broken relations in
Sec. 8.4):
    eq. (8.5.24)-(8.5.31)            Stone C/D/E/F formulas (e.g. eq.8.5.26
                                        exponent 1/3, eq.8.5.27 phase (a-1)/2)

Usage:
    python3 special_values_8_5.py [--n N] [--seed S]
"""
from __future__ import annotations

import argparse
import random

from sympy import Rational, Integer, sqrt, factorial as fac, binomial as binom, S
from sympy.physics.wigner import clebsch_gordan as CG

HALF = Rational(1, 2)
TOL = S(10) ** (-18)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _eq(x, y) -> bool:
    d = (x - y).evalf(30)
    return d.is_finite and abs(d) < TOL


def ph(exponent) -> int:
    return int(S.NegativeOne ** exponent)


def rj(lo, hi):
    """random momentum in {lo, lo+1/2, ..., hi} (half-integer lattice)."""
    return Rational(random.randint(int(2 * lo), int(2 * hi)), 2)


def rproj(j):
    """random projection in {-j, -j+1, ..., j}."""
    return -j + random.randint(0, int(2 * j))


def rc(a, b):
    """random third momentum c in {|a-b|, ..., a+b} (unit steps)."""
    return abs(a - b) + random.randint(0, int(2 * min(a, b)))


# ---------------------------------------------------------------------------
# 8.5.1  Special values of the momenta a, b, c
#
# every sampler returns (a, alpha, b, beta, c, gamma, rhs) or None (redraw)
# ---------------------------------------------------------------------------
def s78():                                    # c=0  ->  a=b, beta=-alpha
    a = rj(HALF, 4); al = rproj(a)
    return (a, al, a, -al, S.Zero, S.Zero, ph(a - al) / sqrt(2 * a + 1))


def s79():                                    # b=0
    a = rj(HALF, 4); al = rproj(a)
    return (a, al, S.Zero, S.Zero, a, al, S.One)


def s80():                                    # c=a+b, alpha=a
    a = rj(HALF, 3); b = rj(HALF, 3); al = a; be = rproj(b)
    rhs = sqrt(fac(2 * a) * fac(2 * b) * fac(a + b + al + be) * fac(a + b - al - be)
               / (fac(2 * a + 2 * b) * fac(a + al) * fac(a - al) * fac(b + be) * fac(b - be)))
    return (a, al, b, be, a + b, al + be, rhs)


def s81():                                    # c=a+b, quasi-binomial (i,n int)
    a = rj(HALF, 3); b = rj(HALF, 3)
    i = random.randint(1, int(2 * b) + 1)
    n = i + random.randint(0, int(2 * a))
    if not (1 <= n <= int(2 * a + 2 * b) + 1):
        return None
    al = a - n + i; be = b - i + 1; ga = a + b - n + 1
    if abs(ga) > a + b:
        return None
    rhs = sqrt(binom(2 * a, n - i) * binom(2 * b, i - 1) / binom(2 * a + 2 * b, n - 1))
    return (a, al, b, be, a + b, ga, rhs)


def s82():                                    # c=a+b, quasi-binomial (2nd form)
    a = rj(HALF, 3); b = rj(HALF, 3)
    i = random.randint(1, int(2 * b) + 1)
    n = random.randint(1, int(2 * a) + 1)
    al = a + n - i; be = -b + i - 1; c = a + b; ga = a - b + n - 1
    if abs(al) > a or abs(be) > b or abs(ga) > c:
        return None
    if not (0 <= 2 * a + n - 1 <= int(2 * a + 2 * b)):
        return None
    rhs = sqrt(binom(2 * a, i - n) * binom(2 * b, i - 1) / binom(2 * a + 2 * b, 2 * a + n - 1))
    return (a, al, b, be, c, ga, rhs)


def s83():                                    # c=a+b, stretched
    a = rj(HALF, 3); b = rj(HALF, 3)
    return (a, a, b, b, a + b, a + b, S.One)


def s84():                                    # c=a+b, top/bottom
    a = rj(HALF, 3); b = rj(HALF, 3)
    return (a, a, b, -b, a + b, a - b, sqrt(fac(2 * a) * fac(2 * b) / fac(2 * a + 2 * b)))


def s85():                                    # c=a+b-1
    a = rj(HALF, 3); b = rj(HALF, 3); c = a + b - 1
    if c < abs(a - b) or c < 0:
        return None
    al = rproj(a); be = rproj(b); ga = al + be
    if abs(ga) > c:
        return None
    rhs = 2 * (b * al - a * be) * sqrt(
        (2 * a + 2 * b - 1) * fac(2 * a - 1) * fac(2 * b - 1)
        * fac(a + b + ga - 1) * fac(a + b - ga - 1)
        / (fac(a + al) * fac(a - al) * fac(b + be) * fac(b - be) * fac(2 * a + 2 * b)))
    return (a, al, b, be, c, ga, rhs)


def s86():                                    # c=a+b-1, zero when alpha/beta=a/b
    a = rj(1, 3); al = rproj(a)
    if abs(2 * al) > 2 * a - 1:
        return None
    return (a, al, a, al, 2 * a - 1, 2 * al, S.Zero)


def s87():                                    # c=a+b-1, particular
    a = rj(HALF, 3); b = rj(HALF, 3); c = a + b - 1
    if c < abs(a - b) or c < 0:
        return None
    return (a, a - 1, b, b, c, c, -sqrt(b / (a + b)))


def s88():                                    # c=a+b-1, particular
    a = rj(HALF, 3); b = rj(HALF, 3); c = a + b - 1
    if c < abs(a - b) or c < 0:
        return None
    return (a, a, b, b - 1, c, c, sqrt(a / (a + b)))


def s89():                                    # c=a+b-1, particular
    # factor is (2a+2b-1); an OCR "2a-2b-1" in Chap8.tex has been corrected.
    a = rj(HALF, 3); b = rj(HALF, 3)
    c = a + b - 1
    if c < abs(a - b) or c < 0:
        return None
    return (a, a, b, -b, c, a - b,
            sqrt(fac(2 * a) * fac(2 * b) * (2 * a + 2 * b - 1) / fac(2 * a + 2 * b)))


def s90():                                    # c=a-b (a>=b)
    b = rj(HALF, 3); a = rj(b, 3); c = a - b
    al = rproj(a); be = rproj(b); ga = al + be
    if abs(ga) > c:
        return None
    rhs = ph(b + be) * sqrt(
        fac(a + al) * fac(a - al) * fac(2 * b) * fac(2 * a - 2 * b + 1)
        / (fac(2 * a + 1) * fac(b + be) * fac(b - be)
           * fac(a - b + ga) * fac(a - b - ga)))
    return (a, al, b, be, c, ga, rhs)


def s91():                                    # c=a-b, quasi-binomial
    b = rj(HALF, 3); a = rj(b, 3)
    i = random.randint(1, int(2 * b) + 1)
    lo = int(2 * b) + 1; hi = int(2 * a) + 1
    if lo > hi:
        return None
    n = random.randint(lo, hi)
    if not (0 <= n - i <= int(2 * a)):
        return None
    al = a - n + i; be = b - i + 1; c = a - b; ga = a + b - n + 1
    if abs(ga) > c:
        return None
    rhs = ph(2 * b - i + 1) * sqrt(
        (2 * a - 2 * b + 1) / (2 * a + 1)
        * binom(2 * b, i - 1) * binom(2 * a - 2 * b, 2 * a - n + 1) / binom(2 * a, n - i))
    return (a, al, b, be, c, ga, rhs)


def s92():                                    # c=a-b, quasi-binomial (2nd form)
    b = rj(HALF, 3); a = rj(b, 3)
    i = random.randint(1 - int(2 * b), 1)
    n = random.randint(1, int(2 * a - 2 * b) + 1)
    al = a - n + i; be = -b - i + 1; c = a - b; ga = a - b - n + 1
    if abs(al) > a or abs(be) > b or abs(ga) > c:
        return None
    if not (0 <= n - i <= int(2 * a)):
        return None
    rhs = ph(i + 1) * sqrt((2 * a - 2 * b + 1) / (2 * a + 1)
                           * binom(2 * b, -i + 1) * binom(2 * a - 2 * b, n - 1) / binom(2 * a, n - i))
    return (a, al, b, be, c, ga, rhs)


def s93():                                    # c=a-b, particular
    b = rj(HALF, 3); a = rj(b, 3)
    return (a, a, b, -b, a - b, a - b, sqrt((2 * a - 2 * b + 1) / (2 * a + 1)))


def s94():                                    # c=a-b+1
    a = rj(HALF, 3); b = rj(HALF, 3); c = a - b + 1
    if c < abs(a - b) or c < 0 or 2 * a - 2 * b + 1 < 0:
        return None
    al = rproj(a); be = rproj(b); ga = al + be
    if abs(ga) > c:
        return None
    rhs = ph(b + be + 1) * 2 * (a * be + b * al + be) * sqrt(
        (2 * a - 2 * b + 3) * fac(2 * b - 1) * fac(2 * a - 2 * b + 1)
        * fac(a + al) * fac(a - al)
        / (fac(2 * a + 2) * fac(b + be) * fac(b - be)
           * fac(a - b + ga + 1) * fac(a - b - ga + 1)))
    return (a, al, b, be, c, ga, rhs)


def s95():                                    # c=a-b+1, particular
    a = rj(HALF, 3); b = rj(HALF, 3); c = a - b + 1
    if c < abs(a - b) or c < 0 or 2 * a - 2 * b + 3 < 0:
        return None
    return (a, a, b, -b, c, a - b,
            sqrt((2 * a - 2 * b + 3) * 2 * b / ((2 * a + 2) * (2 * a + 1))))


def s96():                                    # c=a+b-2
    a = rj(1, 3); b = rj(1, 3); c = a + b - 2
    if c < abs(a - b) or c < 0:
        return None
    al = rproj(a); be = rproj(b); ga = al + be
    if abs(ga) > c:
        return None
    pref = sqrt(2 * a * (2 * a - 1) * 2 * b * (2 * b - 1)
                / (2 * (2 * a + 2 * b - 2) * (2 * a + 2 * b - 1)))
    inv = (binom(2 * a, a - al) * binom(2 * b, b - be)
           * binom(2 * a + 2 * b - 4, a + b - ga - 2)) ** (-HALF)
    brace = (binom(2 * a - 2, a - al) * binom(2 * b - 2, b + be)
             - 2 * binom(2 * a - 2, a - al - 1) * binom(2 * b - 2, b + be - 1)
             + binom(2 * a - 2, a - al - 2) * binom(2 * b - 2, b + be - 2))
    return (a, al, b, be, c, ga, pref * inv * brace)


def s97():                                    # c=a+b-2, particular
    a = rj(1, 3); b = rj(1, 3); c = a + b - 2
    if c < abs(a - b) or c < 0:
        return None
    return (a, a, b, -b, c, a - b,
            sqrt(fac(2 * a) * fac(2 * b) * (2 * a + 2 * b - 3) / (2 * fac(2 * a + 2 * b - 1))))


def s98():                                    # c=a-b+2
    a = rj(HALF, 3); b = rj(1, 3); c = a - b + 2
    if c < abs(a - b) or c < 0 or 2 * a - 2 * b + 3 < 0:
        return None
    al = rproj(a); be = rproj(b); ga = al + be
    if abs(ga) > c:
        return None
    pref = sqrt((2 * b - 1) * 2 * b * (2 * a - 2 * b + 5) * (2 * a - 2 * b + 4)
                * (2 * a - 2 * b + 3) / (2 * (2 * a + 1) * (2 * a + 2) * (2 * a + 3)))
    inv = (binom(2 * a, a - al) * binom(2 * b, b - be)
           * binom(2 * a - 2 * b + 4, a - b - ga + 2)) ** (-HALF)
    brace = (binom(2 * b - 2, b + be) * binom(2 * a - 2 * b + 2, a - b - ga)
             - 2 * binom(2 * b - 2, b + be - 1) * binom(2 * a - 2 * b + 2, a - b - ga + 1)
             + binom(2 * b - 2, b + be - 2) * binom(2 * a - 2 * b + 2, a - b - ga + 2))
    return (a, al, b, be, c, ga, pref * ph(b + be) * inv * brace)


def s99():                                    # c=a-b+2, particular
    a = rj(HALF, 3); b = rj(1, 3); c = a - b + 2
    if c < abs(a - b) or c < 0 or 2 * a - 2 * b + 5 < 0:
        return None
    return (a, a, b, -b, c, a - b,
            sqrt((2 * a - 2 * b + 5) * 2 * b * (2 * b - 1)
                 / ((2 * a + 1) * (2 * a + 2) * (2 * a + 3))))


def s100():                                   # a=b, alpha=beta ; gamma=2 alpha
    a = rj(HALF, 3); al = rproj(a); ga = 2 * al
    lo = int(abs(2 * al)); hi = int(2 * a)
    if lo > hi:
        return None
    c = random.randint(lo, hi)
    g2 = int(2 * a) + c
    if g2 % 2 == 1:                           # 2a+c = 2g+1  ->  0
        rhs = S.Zero
    else:                                     # 2a+c = 2g
        g = Integer(g2 // 2)
        cp = (Integer(c) + ga); cm = (Integer(c) - ga)
        if not (cp / 2).is_integer or not (cm / 2).is_integer or g - c < 0:
            return None
        rhs = (ph(g - c) * sqrt(2 * c + 1) * fac(g)
               / (fac(cp / 2) * fac(cm / 2) * fac(g - c))
               * sqrt(fac(cp) * fac(cm) * fac(2 * g - 2 * c) / fac(2 * g + 1)))
    return (a, al, a, al, Integer(c), ga, rhs)


# ---------------------------------------------------------------------------
# 8.5.2  Special values of the momentum projections
# ---------------------------------------------------------------------------
def s109():                                   # alpha=beta=gamma=0  (a,b,c integer)
    a = Integer(random.randint(0, 3)); b = Integer(random.randint(0, 3))
    c = Integer(random.randint(int(abs(a - b)), int(a + b)))
    if (a + b + c) % 2 == 1:
        rhs = S.Zero
    else:
        g = Rational(a + b + c, 2)
        rhs = (ph(g - c) * sqrt(2 * c + 1) * fac(g) / (fac(g - a) * fac(g - b) * fac(g - c))
               * sqrt(fac(2 * g - 2 * a) * fac(2 * g - 2 * b) * fac(2 * g - 2 * c) / fac(2 * g + 1)))
    return (a, S.Zero, b, S.Zero, c, S.Zero, rhs)


def s110():                                   # (a0 b0 | a+b 0)
    a = Integer(random.randint(0, 3)); b = Integer(random.randint(0, 3))
    return (a, S.Zero, b, S.Zero, a + b, S.Zero,
            fac(a + b) / (fac(a) * fac(b)) * sqrt(fac(2 * a) * fac(2 * b) / fac(2 * a + 2 * b)))


def s111():                                   # (a0 b0 | a-b 0)
    b = Integer(random.randint(0, 3)); a = Integer(random.randint(int(b), 3))
    return (a, S.Zero, b, S.Zero, a - b, S.Zero,
            ph(b) * fac(a) / (fac(b) * fac(a - b))
            * sqrt(fac(2 * b) * fac(2 * a - 2 * b + 1) / fac(2 * a + 1)))


def s112():                                   # gamma=c
    a = rj(HALF, 3); b = rj(HALF, 3); c = rc(a, b)
    al = rproj(a); be = c - al
    if abs(be) > b:
        return None
    rhs = ph(a - al) * sqrt(
        fac(2 * c + 1) * fac(a + b - c) * fac(a + al) * fac(b + be)
        / (fac(a + b + c + 1) * fac(a - b + c) * fac(-a + b + c) * fac(a - al) * fac(b - be)))
    return (a, al, b, be, c, c, rhs)


def s113():                                   # alpha=a
    a = rj(HALF, 3); b = rj(HALF, 3); c = rc(a, b)
    be = rproj(b); ga = a + be
    if abs(ga) > c:
        return None
    rhs = sqrt(
        (2 * c + 1) * fac(2 * a) * fac(-a + b + c) * fac(b - be) * fac(c + ga)
        / (fac(a + b + c + 1) * fac(a - b + c) * fac(a + b - c) * fac(b + be) * fac(c - ga)))
    return (a, a, b, be, c, ga, rhs)


def s114():                                   # (aa bb | cc) = delta_{a+b,c}
    a = rj(HALF, 3); b = rj(HALF, 3); c = rc(a, b)
    return (a, a, b, b, c, c, S.One if c == a + b else S.Zero)


def s115():                                   # (aa b-b | cc)
    a = rj(HALF, 3); b = rj(HALF, 3); c = rc(a, b)
    rhs = sqrt((2 * c + 1) / (2 * a + 1)) if c == a - b else S.Zero
    return (a, a, b, -b, c, c, rhs)


def s116():                                   # (aa b,c-a | cc)
    a = rj(HALF, 3); b = rj(HALF, 3); c = rc(a, b); be = c - a
    if abs(be) > b:
        return None
    return (a, a, b, be, c, c,
            sqrt(fac(2 * a) * fac(2 * c + 1) / (fac(a + b + c + 1) * fac(a - b + c))))


def s117():                                   # (a,a-1, b,c-a+1 | cc)
    a = rj(HALF, 3); b = rj(HALF, 3); c = rc(a, b); be = c - a + 1
    if abs(be) > b:
        return None
    return (a, a - 1, b, be, c, c,
            -sqrt(fac(2 * a - 1) * fac(2 * c + 1) * (a + b - c) * (-a + b + c + 1)
                  / (fac(a + b + c + 1) * fac(a - b + c))))


def s118():                                   # (aa b,c-a-1 | c,c-1)
    a = rj(HALF, 3); b = rj(HALF, 3); c = rc(a, b)
    if c < 1:
        return None
    be = c - a - 1
    if abs(be) > b:
        return None
    return (a, a, b, be, c, c - 1,
            sqrt(fac(2 * a) * fac(2 * c + 1) * (a + b - c + 1) * (-a + b + c)
                 / (fac(a + b + c + 1) * fac(a - b + c) * 2 * c)))


def s119():                                   # (cc b0 | cc)
    c = rj(HALF, 3); b = Integer(random.randint(0, int(2 * c)))
    return (c, c, b, S.Zero, c, c,
            fac(2 * c) * sqrt((2 * c + 1) / (fac(2 * c - b) * fac(2 * c + b + 1))))


def s120():                                   # (c,c-b, bb | cc)
    # numerator is the factorial (2c+1)!; an OCR "(2c+1)" in Chap8.tex has been
    # corrected.
    c = rj(HALF, 3); b = Integer(random.randint(0, int(2 * c)))
    return (c, c - b, b, b, c, c,
            ph(b) * sqrt(fac(2 * c + 1) * fac(2 * b) / (fac(2 * c + b + 1) * fac(b))))


def s121():                                   # gamma=c-1
    a = rj(HALF, 3); b = rj(HALF, 3); c = rc(a, b)
    if abs(c - 1) > c:
        return None
    al = rproj(a); be = (c - 1) - al
    if abs(be) > b:
        return None
    brace = (b - be) * (b + be + 1) - (a - al) * (a + al + 1)
    rhs = ph(a - al) * brace * sqrt(
        (2 * c + 1) * fac(2 * c - 1) * fac(a + b - c) * fac(a + al) * fac(b + be)
        / (fac(a + b + c + 1) * fac(a - b + c) * fac(-a + b + c) * fac(a - al) * fac(b - be)))
    return (a, al, b, be, c, c - 1, rhs)


def s122():                                   # alpha=a-1
    a = rj(HALF, 3); b = rj(HALF, 3); c = rc(a, b)
    be = rproj(b); ga = (a - 1) + be
    if abs(ga) > c:
        return None
    brace = (c - ga) * (c + ga + 1) - (b + be) * (b - be + 1)
    rhs = brace * sqrt(
        (2 * c + 1) * fac(2 * a - 1) * fac(-a + b + c) * fac(b - be) * fac(c + ga)
        / (fac(a + b + c + 1) * fac(a - b + c) * fac(a + b - c) * fac(b + be) * fac(c - ga)))
    return (a, a - 1, b, be, c, ga, rhs)


def s123():                                   # (a al, a al | c,c-1) = 0
    # Two identical (a,alpha) coupled to (c,c-1) vanishes identically: exchange
    # gives a factor (-1)^{2a-c}, and 2a-c is forced odd here (c-1 = 2 alpha).
    a = rj(HALF, 3); al = rproj(a); c = rc(a, a)
    if abs(c - 1) > c:
        return None
    return (a, al, a, al, c, c - 1, S.Zero)


def s124():                                   # (a,a-1, c,-gamma | c,gamma) = 0
    a = Integer(random.choice([1, 3])); ga = Rational(a - 1, 2)
    cc = Integer(random.randint(int(a), 3))
    return (a, a - 1, cc, -ga, cc, ga, S.Zero)


def s125():                                   # (a,a-1, b,beta | c,c-1)
    a = rj(HALF, 3); b = rj(HALF, 3); c = rc(a, b)
    if abs(c - 1) > c:
        return None
    be = c - a
    if abs(be) > b:
        return None
    brace = a * (a + 1) - b * (b + 1) + c * (c + 1) - 2 * a * c
    rhs = brace * sqrt(
        (2 * c + 1) * fac(2 * c - 1) * fac(2 * a - 1) / (fac(a + b + c + 1) * fac(a - b + c)))
    return (a, a - 1, b, be, c, c - 1, rhs)


# ---------------------------------------------------------------------------
# registry
# ---------------------------------------------------------------------------
SEC_851 = [
    ("eq 8.5.1   c=0", s78),
    ("eq 8.5.2   b=0", s79),
    ("eq 8.5.3   c=a+b, alpha=a", s80),
    ("eq 8.5.4   c=a+b (quasi-binom)", s81),
    ("eq 8.5.5   c=a+b (quasi-binom 2)", s82),
    ("eq 8.5.6   c=a+b, stretched", s83),
    ("eq 8.5.7   c=a+b, (aa,b-b)", s84),
    ("eq 8.5.8   c=a+b-1", s85),
    ("eq 8.5.9   c=a+b-1, zero rule", s86),
    ("eq 8.5.10   c=a+b-1, particular", s87),
    ("eq 8.5.11   c=a+b-1, particular", s88),
    ("eq 8.5.12   c=a+b-1, particular", s89),
    ("eq 8.5.13   c=a-b", s90),
    ("eq 8.5.14   c=a-b (quasi-binom)", s91),
    ("eq 8.5.15   c=a-b (quasi-binom 2)", s92),
    ("eq 8.5.16   c=a-b, particular", s93),
    ("eq 8.5.17   c=a-b+1", s94),
    ("eq 8.5.18   c=a-b+1, particular", s95),
    ("eq 8.5.19   c=a+b-2", s96),
    ("eq 8.5.20   c=a+b-2, particular", s97),
    ("eq 8.5.21   c=a-b+2", s98),
    ("eq 8.5.22   c=a-b+2, particular", s99),
    ("eq 8.5.23  a=b, alpha=beta", s100),
]

SEC_852 = [
    ("eq 8.5.32  alpha=beta=gamma=0", s109),
    ("eq 8.5.33  (a0 b0|a+b 0)", s110),
    ("eq 8.5.34  (a0 b0|a-b 0)", s111),
    ("eq 8.5.35  gamma=c", s112),
    ("eq 8.5.36  alpha=a", s113),
    ("eq 8.5.37  (aa bb|cc)", s114),
    ("eq 8.5.38  (aa b-b|cc)", s115),
    ("eq 8.5.39  (aa b,c-a|cc)", s116),
    ("eq 8.5.40  (a,a-1,...|cc)", s117),
    ("eq 8.5.41  (aa...|c,c-1)", s118),
    ("eq 8.5.42  (cc b0|cc)", s119),
    ("eq 8.5.43  (c,c-b,bb|cc)", s120),
    ("eq 8.5.44  gamma=c-1", s121),
    ("eq 8.5.45  alpha=a-1", s122),
    ("eq 8.5.46  zero rule", s123),
    ("eq 8.5.47  zero rule", s124),
    ("eq 8.5.48  (a,a-1,...|c,c-1)", s125),
]


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def run(n, seed):
    random.seed(seed)
    print(f"Section 8.5 special-value check -- seed={seed}, up to {n} random "
          f"instances per formula\n")

    all_ok = True
    for title, table in (("8.5.1  Special values of momenta a,b,c", SEC_851),
                         ("8.5.2  Special values of projections", SEC_852)):
        print(f"=== {title} ===")
        for label, sampler in table:
            got, nonzero, bad = 0, 0, None
            for _ in range(n * 400):
                if got >= n:
                    break
                draw = sampler()
                if draw is None:
                    continue
                a, al, b, be, c, ga, rhs = draw
                lhs = CG(a, b, c, al, be, ga)
                got += 1
                if not _eq(lhs, rhs):
                    if bad is None:
                        bad = ((a, al, b, be, c, ga), lhs, rhs)
                elif lhs != 0:
                    nonzero += 1
            if got == 0:
                print(f"  [SKIP] {label:34s} no valid draws")
                continue
            ok = bad is None
            all_ok &= ok
            mark = "OK  " if ok else "FAIL"
            print(f"  [{mark}] {label:34s} {got} instances ({nonzero} non-zero)")
            if bad is not None:
                cfg, lhs, rhs = bad
                print(f"         counterexample (a,al,b,be,c,ga)={cfg}")
                print(f"           lhs(CG)={lhs}   rhs(book)={rhs}")
        print()

    print("Not tested (auxiliary-function formulas with OCR artifacts -- flag for repair):")
    print("    eq. (8.5.24)-(8.5.31)     Stone C/D/E/F formulas")
    print()
    print("ALL FORMULAS HOLD" if all_ok else "SOME FORMULAS FAILED -- see above")
    return all_ok


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Check Section 8.5 special-value formulas.")
    p.add_argument("--n", type=int, default=6, help="instances per formula")
    p.add_argument("--seed", type=int, default=20260805, help="RNG seed")
    args = p.parse_args()
    raise SystemExit(0 if run(args.n, args.seed) else 1)
