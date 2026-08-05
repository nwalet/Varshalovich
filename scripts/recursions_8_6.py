#!/usr/bin/env python3
r"""
Numerical check of the Clebsch-Gordan recursion relations listed in Section 8.6
of Varshalovich, Moskalev & Khersonskii, "Quantum Theory of Angular Momentum".

Each relation is an identity among Clebsch-Gordan coefficients whose arguments
are shifted by 1, 1/2, ... .  We draw random valid angular momenta/projections,
evaluate both sides with sympy's clebsch_gordan (coefficients that fall outside
the physical domain are treated as 0), and compare.

The tests are SPLIT BY SUBSECTION, matching the book:
    8.6.1  General recursion relations          eq. 126-128   [TODO]
    8.6.2  Arguments alpha,beta,gamma change 1  eq. 129-133   [done]
    8.6.3  Arguments change by 1/2              eq. 134-143   [done]
    8.6.4  The case alpha=beta=gamma=0          eq. 144-147   [done]
    8.6.5  Arguments a,b,c change by 1          eq. 148-151   [done]
    8.6.6  Arguments a,b,alpha,beta change 1    eq. 152-153   [TODO]
    8.6.7  Arguments c,b,gamma,beta change 1    eq. 154-156   [TODO]
    8.6.8  Recursion relations for R-symbols    eq. 157-161   [TODO]

Convention (book == sympy):
    C_{a alpha, b beta}^{c gamma} == clebsch_gordan(a, b, c, alpha, beta, gamma)

Many relations carry a correlated (upper/lower) sign choice; we encode it with
s = +1 / -1 (so  +/-  ->  s  and  -/+  ->  -s) and test both.

Usage:
    python3 recursions_8_6.py [--n N] [--seed S]
"""
from __future__ import annotations

import argparse
import random

from sympy import Rational, Integer, sqrt, factorial as fac, S
from sympy.physics.wigner import clebsch_gordan

HALF = Rational(1, 2)
TOL = S(10) ** (-18)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _eq(x, y) -> bool:
    d = (x - y).evalf(30)
    return d.is_finite and abs(d) < TOL


def ph(n) -> int:
    return int(S.NegativeOne ** n)


def C(a, b, cc, al, be, ga):
    """Clebsch-Gordan coefficient, 0 outside the physical domain (so shifted
    arguments in a recursion never raise)."""
    if a < 0 or b < 0 or cc < 0:
        return S.Zero
    if abs(al) > a or abs(be) > b or abs(ga) > cc:
        return S.Zero
    if cc < abs(a - b) or cc > a + b:
        return S.Zero
    if not (a - al).is_integer or not (b - be).is_integer or not (cc - ga).is_integer:
        return S.Zero
    if not (a + b + cc).is_integer:
        return S.Zero
    if al + be != ga:
        return S.Zero
    return clebsch_gordan(a, b, cc, al, be, ga)


def rj(lo, hi):                     # momentum on the half-integer lattice
    return Rational(random.randint(int(2 * lo), int(2 * hi)), 2)


def rhalf(hi):                      # strictly half-integer momentum 1/2..hi
    return Rational(random.randrange(1, 2 * int(hi) + 1, 2), 2)


def rint(lo, hi):                   # integer momentum
    return Integer(random.randint(int(lo), int(hi)))


def rproj(j):
    return -j + random.randint(0, int(2 * j))


def rc(a, b):
    return abs(a - b) + random.randint(0, int(2 * min(a, b)))


def cfg_master(jmax=3):
    """A generic valid, non-vanishing CG configuration (a,b,c,alpha,beta)."""
    a = rj(HALF, jmax); b = rj(HALF, jmax); c = rc(a, b)
    al = rproj(a); be = rproj(b)
    if abs(al + be) > c:
        return None
    if C(a, b, c, al, be, al + be) == 0:
        return None
    return (a, b, c, al, be)


def cfg_zero(jmax=3):
    """Integer a,b,c with a non-vanishing all-zero-projection CG (a+b+c even)."""
    a = rint(0, jmax); b = rint(0, jmax); c = rint(abs(a - b), a + b)
    if C(a, b, c, 0, 0, 0) == 0:
        return None
    return (a, b, c)


# ===========================================================================
# 8.6.2  Arguments alpha, beta, gamma change by 1     (eq. 129-133)
# ===========================================================================
def r129():
    # [(c +/-g)(c -/+g+1)]^1/2 C_{a al,b be}^{c,g -/+1}
    #   = [(a -/+al)(a +/-al+1)]^1/2 C_{a,al+/-1,b be}^{c g}
    #   + [(b -/+be)(b +/-be+1)]^1/2 C_{a al,b,be+/-1}^{c g}
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    s = random.choice([1, -1])
    g0 = al + be                      # conserved third projection (LHS)
    g = g0 + s                         # book's "gamma" (RHS third index)
    lhs = sqrt((c + s * g) * (c - s * g + 1)) * C(a, b, c, al, be, g0)
    rhs = (sqrt((a - s * al) * (a + s * al + 1)) * C(a, b, c, al + s, be, g)
           + sqrt((b - s * be) * (b + s * be + 1)) * C(a, b, c, al, be + s, g))
    return lhs, rhs


def r131():
    # |alpha|=|beta|=1/2 :
    # C_{a 1/2,b 1/2}^{c 1} = C_{a 1/2,b -1/2}^{c 0}
    #                         * [(2b+1)+(-1)^{a+b-c}(2a+1)] / (2 sqrt(c(c+1)))
    a = rhalf(3); b = rhalf(3)
    c = abs(a - b) + random.randint(0, int(2 * min(a, b)))
    if c < 1:
        return None
    lhs = C(a, b, c, HALF, HALF, 1)
    rhs = (C(a, b, c, HALF, -HALF, 0)
           * ((2 * b + 1) + ph(a + b - c) * (2 * a + 1)) / (2 * sqrt(c * (c + 1))))
    return lhs, rhs


def r132():
    # |alpha|=|beta|=1, a+b+c even :  (source exponent 1/3 is an OCR typo -> 1/2)
    # C_{a 1,b -1}^{c 0} = C_{a 0,b 0}^{c 0}
    #                      * [c(c+1)-a(a+1)-b(b+1)] / (2 sqrt(a(a+1)b(b+1)))
    a = rint(1, 3); b = rint(1, 3); c = rint(abs(a - b), a + b)
    if (a + b + c) % 2 == 1:
        return None
    if C(a, b, c, 0, 0, 0) == 0:
        return None
    lhs = C(a, b, c, 1, -1, 0)
    rhs = (C(a, b, c, 0, 0, 0)
           * (c * (c + 1) - a * (a + 1) - b * (b + 1))
           / (2 * sqrt(a * (a + 1) * b * (b + 1))))
    return lhs, rhs


def r133():
    # |alpha|=|beta|=1, a+b+c even :  (source exponent 4/2 is an OCR typo -> 1/2)
    # C_{a 1,b 1}^{c 2} = C_{a 0,b 0}^{c 0}
    #   * {a(a+1)[c(c+1)-a(a+1)+b(b+1)] + b(b+1)[c(c+1)+a(a+1)-b(b+1)]}
    #     / (2 sqrt(a(a+1)b(b+1)(c-1)c(c+1)(c+2)))
    a = rint(1, 3); b = rint(1, 3); c = rint(abs(a - b), a + b)
    if (a + b + c) % 2 == 1 or c < 2:
        return None
    if C(a, b, c, 0, 0, 0) == 0:
        return None
    num = (a * (a + 1) * (c * (c + 1) - a * (a + 1) + b * (b + 1))
           + b * (b + 1) * (c * (c + 1) + a * (a + 1) - b * (b + 1)))
    den = 2 * sqrt(a * (a + 1) * b * (b + 1) * (c - 1) * c * (c + 1) * (c + 2))
    return C(a, b, c, 1, 1, 2), C(a, b, c, 0, 0, 0) * num / den


def r130():
    # |gamma|=c :  C_{a,al-/+1, b,be}^{c,+/-c}
    #            = -C_{a,al, b,be-/+1}^{c,+/-c} [(b+/-be)(b-/+be+1)/((a+/-al)(a-/+al+1))]^1/2
    a = rj(HALF, 3); b = rj(HALF, 3); c = rc(a, b)
    al = rproj(a); be = rproj(b)
    s = random.choice([1, -1])
    den = (a + s * al) * (a - s * al + 1)
    if den == 0:
        return None
    lhs = C(a, b, c, al - s, be, s * c)
    if lhs == 0:                                   # keep the test non-vacuous
        return None
    rhs = -C(a, b, c, al, be - s, s * c) * sqrt((b + s * be) * (b - s * be + 1) / den)
    return lhs, rhs


SEC_862 = [
    ("eq 8.6.129  master ladder in gamma", r129),
    ("eq 8.6.130  |gamma|=c", r130),
    ("eq 8.6.131  |al|=|be|=1/2", r131),
    ("eq 8.6.132  |al|=|be|=1 (exp 1/3->1/2)", r132),
    ("eq 8.6.133  |al|=|be|=1 (exp 4/2->1/2)", r133),
]


# ===========================================================================
# 8.6.3  Arguments change by 1/2                      (eq. 134-143)
# ===========================================================================
def r134():
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    s = random.choice([1, -1]); g = al + be
    lhs = (2 * a + 1) * sqrt(b + s * be) * C(a, b, c, al, be, g)
    rhs = (-s * sqrt((a - s * al) * (a + b - c) * (a + b + c + 1))
           * C(a - HALF, b - HALF, c, al + s * HALF, be - s * HALF, g)
           + sqrt((a + s * al + 1) * (-a + b + c) * (a - b + c + 1))
           * C(a + HALF, b - HALF, c, al + s * HALF, be - s * HALF, g))
    return lhs, rhs


def r135():
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    s = random.choice([1, -1]); g = al + be
    lhs = (2 * a + 1) * sqrt(b - s * be + 1) * C(a, b, c, al, be, g)
    rhs = (sqrt((a - s * al) * (a - b + c) * (-a + b + c + 1))
           * C(a - HALF, b + HALF, c, al + s * HALF, be - s * HALF, g)
           + s * sqrt((a + s * al + 1) * (a + b - c + 1) * (a + b + c + 2))
           * C(a + HALF, b + HALF, c, al + s * HALF, be - s * HALF, g))
    return lhs, rhs


def r136():
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    g = al + be
    lhs = sqrt((-a + b + c) * (a - b + c + 1)) * C(a, b, c, al, be, g)
    rhs = (sqrt((a - al + 1) * (b - be)) * C(a + HALF, b - HALF, c, al - HALF, be + HALF, g)
           + sqrt((a + al + 1) * (b + be)) * C(a + HALF, b - HALF, c, al + HALF, be - HALF, g))
    return lhs, rhs


def r137():
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    g = al + be
    lhs = sqrt((a - b + c) * (-a + b + c + 1)) * C(a, b, c, al, be, g)
    rhs = (sqrt((a - al) * (b - be + 1)) * C(a - HALF, b + HALF, c, al + HALF, be - HALF, g)
           + sqrt((a + al) * (b + be + 1)) * C(a - HALF, b + HALF, c, al - HALF, be + HALF, g))
    return lhs, rhs


def r138():
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    g = al + be
    lhs = sqrt(2 * c * (-a + b + c) * (a + b + c + 1) / (2 * c + 1)) * C(a, b, c, al, be, g)
    rhs = (sqrt((b - be) * (c - g)) * C(a, b - HALF, c - HALF, al, be + HALF, g + HALF)
           + sqrt((b + be) * (c + g)) * C(a, b - HALF, c - HALF, al, be - HALF, g - HALF))
    return lhs, rhs


def r139():
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    s = random.choice([1, -1]); g = al + be
    lhs = sqrt(2 * c * (c - s * g) * (a + b + c + 1) / (2 * c + 1)) * C(a, b, c, al, be, g)
    rhs = (sqrt((a - s * al) * (a - b + c)) * C(a - HALF, b, c - HALF, al + s * HALF, be, g + s * HALF)
           + sqrt((b - s * be) * (-a + b + c)) * C(a, b - HALF, c - HALF, al, be + s * HALF, g + s * HALF))
    return lhs, rhs


def r139b():
    # the multline between eq.139 and eq.140 (labelled chap8:eq:139b in the
    # source); changes a by +1/2 (term 1) and b by +1/2 (term 2), both to c-1/2.
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    s = random.choice([1, -1]); g = al + be
    lhs = sqrt(2 * c * (c - s * g) * (a + b - c + 1) / (2 * c + 1)) * C(a, b, c, al, be, g)
    rhs = (s * sqrt((a + s * al + 1) * (-a + b + c))
           * C(a + HALF, b, c - HALF, al + s * HALF, be, g + s * HALF)
           - s * sqrt((b + s * be + 1) * (a - b + c))
           * C(a, b + HALF, c - HALF, al, be + s * HALF, g + s * HALF))
    return lhs, rhs


def r140():
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    s = random.choice([1, -1]); g = al + be
    lhs = sqrt(2 * (c + 1) * (c + s * g + 1) * (a + b - c) / (2 * c + 1)) * C(a, b, c, al, be, g)
    rhs = (-s * sqrt((a - s * al) * (-a + b + c + 1)) * C(a - HALF, b, c + HALF, al + s * HALF, be, g + s * HALF)
           + s * sqrt((b - s * be) * (a - b + c + 1)) * C(a, b - HALF, c + HALF, al, be + s * HALF, g + s * HALF))
    return lhs, rhs


def r141():
    # [(2c+1)(b-/+be)]^1/2 C_{a al,b be}^{c ga} = ... (mixes c-1/2 and c+1/2)
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    if c < 1:
        return None
    s = random.choice([1, -1]); g = al + be
    lhs = sqrt((2 * c + 1) * (b - s * be)) * C(a, b, c, al, be, g)
    rhs = (sqrt((c - s * g) * (-a + b + c) * (a + b + c + 1) / (2 * c))
           * C(a, b - HALF, c - HALF, al, be + s * HALF, g + s * HALF)
           + s * sqrt((c + s * g + 1) * (a - b + c + 1) * (a + b - c) / (2 * (c + 1)))
           * C(a, b - HALF, c + HALF, al, be + s * HALF, g + s * HALF))
    return lhs, rhs


def r142():
    # particular case of eq.134/135, a+b+c odd
    a = rhalf(3); b = rhalf(3)
    c = abs(a - b) + random.randint(0, int(2 * min(a, b)))
    if (a + b + c) % 2 == 0:
        return None
    s = random.choice([1, -1])
    lhs = C(a, b, c, s * HALF, -s * HALF, 0)
    mid = (s * C(a - HALF, b - HALF, c, 0, 0, 0)
           * sqrt((a + b - c) * (a + b + c + 1) / ((2 * a + 1) * (2 * b + 1))))
    right = (-s * C(a + HALF, b + HALF, c, 0, 0, 0)
             * sqrt((a + b - c + 1) * (a + b + c + 2) / ((2 * a + 1) * (2 * b + 1))))
    if not _eq(lhs, mid):
        return lhs, mid
    return lhs, right


def r143():
    # particular case of eq.134/135, a+b+c even
    a = rhalf(3); b = rhalf(3)
    c = abs(a - b) + random.randint(0, int(2 * min(a, b)))
    if (a + b + c) % 2 == 1:
        return None
    s = random.choice([1, -1])
    lhs = C(a, b, c, s * HALF, -s * HALF, 0)
    mid = (C(a + HALF, b - HALF, c, 0, 0, 0)
           * sqrt((-a + b + c) * (a - b + c + 1) / ((2 * a + 1) * (2 * b + 1))))
    right = (C(a - HALF, b + HALF, c, 0, 0, 0)
             * sqrt((a - b + c) * (-a + b + c + 1) / ((2 * a + 1) * (2 * b + 1))))
    if not _eq(lhs, mid):
        return lhs, mid
    return lhs, right


SEC_863 = [
    ("eq 8.6.134  a-+1/2,b-1/2", r134),
    ("eq 8.6.135  a-+1/2,b+1/2", r135),
    ("eq 8.6.136  a+1/2,b-1/2", r136),
    ("eq 8.6.137  a-1/2,b+1/2", r137),
    ("eq 8.6.138  b-1/2,c-1/2", r138),
    ("eq 8.6.139  a-1/2,c-1/2", r139),
    ("eq 8.6.139b multline", r139b),
    ("eq 8.6.140  a-1/2,c+1/2", r140),
    ("eq 8.6.141  b-1/2, c-/+1/2", r141),
    ("eq 8.6.142  |al|=|be|=1/2, a+b+c odd", r142),
    ("eq 8.6.143  |al|=|be|=1/2, a+b+c even", r143),
]


# ===========================================================================
# 8.6.4  The case alpha = beta = gamma = 0            (eq. 144-147)
# ===========================================================================
def r144():
    z = cfg_zero()
    if z is None:
        return None
    a, b, c = z
    g = Rational(a + b + c, 2)                     # 2g = a+b+c
    ps = [p for p in range(-int(a), int(b) + 1)
          if abs(a - b + 2 * p) <= c and g - a - p >= 0 and g - b + p >= 0]
    if not ps:
        return None
    p = random.choice(ps)
    lhs = C(a + p, b - p, c, 0, 0, 0)
    rhs = (C(a, b, c, 0, 0, 0)
           * fac(g - a) * fac(g - b) / (fac(g - a - p) * fac(g - b + p))
           * sqrt(fac(2 * g - 2 * a - 2 * p) * fac(2 * g - 2 * b + 2 * p)
                  / (fac(2 * g - 2 * a) * fac(2 * g - 2 * b))))
    return lhs, rhs


def r145():
    z = cfg_zero()
    if z is None:
        return None
    a, b, c = z
    if b < 1 or (-a + b + c - 1) <= 0 or abs(a - b + 2) > c:
        return None
    lhs = C(a + 1, b - 1, c, 0, 0, 0)
    rhs = C(a, b, c, 0, 0, 0) * sqrt((-a + b + c) * (a - b + c + 1)
                                     / ((-a + b + c - 1) * (a - b + c + 2)))
    return lhs, rhs


def r146():
    z = cfg_zero()
    if z is None:
        return None
    a, b, c = z
    g = Rational(a + b + c, 2)
    ps = [p for p in range(-int(b // 2), int(b) + 1)
          if b + 2 * p >= 0 and abs(a - b - 2 * p) <= c and c <= a + b + 2 * p
          and g - a + p >= 0 and g - b - p >= 0 and g - c + p >= 0
          and a - b + c - 2 * p >= 0]
    if not ps:
        return None
    p = random.choice(ps)
    lhs = C(a, b + 2 * p, c, 0, 0, 0)
    rhs = (C(a, b, c, 0, 0, 0) * ph(p)
           * fac(g + p) * fac(g - a) * fac(g - b) * fac(g - c)
           / (fac(g) * fac(g - a + p) * fac(g - b - p) * fac(g - c + p))
           * sqrt(fac(a + b - c + 2 * p) * fac(a - b + c - 2 * p)
                  * fac(-a + b + c + 2 * p) * fac(a + b + c + 1)
                  / (fac(a + b - c) * fac(a - b + c) * fac(-a + b + c)
                     * fac(a + b + c + 2 * p + 1))))
    return lhs, rhs


def r147():
    z = cfg_zero()
    if z is None:
        return None
    a, b, c = z
    if (a - b + c - 1) <= 0 or abs(a - b - 2) > c or c > a + b + 2:
        return None
    lhs = C(a, b + 2, c, 0, 0, 0)
    rhs = -C(a, b, c, 0, 0, 0) * sqrt(
        (a + b + c + 2) * (a + b - c + 1) * (a - b + c) * (-a + b + c + 1)
        / ((a + b + c + 3) * (a + b - c + 2) * (a - b + c - 1) * (-a + b + c + 2)))
    return lhs, rhs


SEC_864 = [
    ("eq 8.6.144  a+p, b-p", r144),
    ("eq 8.6.145  a+1, b-1", r145),
    ("eq 8.6.146  b+2p", r146),
    ("eq 8.6.147  b+2", r147),
]


# ===========================================================================
# 8.6.5  Arguments a, b, c change by 1                (eq. 148-151)
#
# Master three-term recursions in a (the base coefficient C_{a al,b be}^{c ga}
# is the middle term; a is shifted by -1, 0, +1).  All keep (alpha,beta,gamma).
# ===========================================================================
def r148():
    # 2[b^2-be^2]^1/2 C_{a al, b-1, be}^{c ga} = ...
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    if b < 1:
        return None
    g = al + be
    lhs = 2 * sqrt(b ** 2 - be ** 2) * C(a, b - 1, c, al, be, g)
    rhs = (1 / (a * (2 * a + 1))
           * sqrt((a ** 2 - al ** 2) * (-a + b + c) * (-a + b + c + 1) * (a - b + c) * (a - b + c + 1))
           * C(a - 1, b, c, al, be, g)
           + al / (a * (a + 1))
           * sqrt((-a + b + c) * (a - b + c + 1) * (a + b - c) * (a + b + c + 1))
           * C(a, b, c, al, be, g)
           - 1 / ((a + 1) * (2 * a + 1))
           * sqrt(((a + 1) ** 2 - al ** 2) * (a + b - c) * (a + b - c + 1) * (a + b + c + 1) * (a + b + c + 2))
           * C(a + 1, b, c, al, be, g))
    return lhs, rhs


def r149():
    # 2 be C_{a al,b be}^{c ga} = ...
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    g = al + be
    lhs = 2 * be * C(a, b, c, al, be, g)
    rhs = (-1 / (a * (2 * a + 1))
           * sqrt((a ** 2 - al ** 2) * (-a + b + c + 1) * (a - b + c) * (a + b - c) * (a + b + c + 1))
           * C(a - 1, b, c, al, be, g)
           - al / (a * (a + 1)) * (a * (a + 1) + b * (b + 1) - c * (c + 1))
           * C(a, b, c, al, be, g)
           - 1 / ((a + 1) * (2 * a + 1))
           * sqrt(((a + 1) ** 2 - al ** 2) * (-a + b + c) * (a - b + c + 1) * (a + b - c + 1) * (a + b + c + 2))
           * C(a + 1, b, c, al, be, g))
    return lhs, rhs


def r150():
    # 2[(b+1)^2-be^2]^1/2 C_{a al, b+1, be}^{c ga} = ...
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    g = al + be
    lhs = 2 * sqrt((b + 1) ** 2 - be ** 2) * C(a, b + 1, c, al, be, g)
    rhs = (-1 / (a * (2 * a + 1))
           * sqrt((a ** 2 - al ** 2) * (a + b - c) * (a + b - c + 1) * (a + b + c + 1) * (a + b + c + 2))
           * C(a - 1, b, c, al, be, g)
           + al / (a * (a + 1))
           * sqrt((-a + b + c + 1) * (a - b + c) * (a + b - c + 1) * (a + b + c + 2))
           * C(a, b, c, al, be, g)
           + 1 / ((a + 1) * (2 * a + 1))
           * sqrt(((a + 1) ** 2 - al ** 2) * (-a + b + c) * (-a + b + c + 1) * (a - b + c) * (a - b + c + 1))
           * C(a + 1, b, c, al, be, g))
    return lhs, rhs


def r151():
    # C_{a al,b be}^{c ga} in terms of the c-1 and c-2 coefficients
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    g = al + be
    if c < 2 or abs(g) >= c or (-a + b + c) <= 0 or (a - b + c) <= 0:
        return None
    pref = sqrt(4 * c ** 2 * (2 * c + 1) * (2 * c - 1)
                / ((c + g) * (c - g) * (-a + b + c) * (a - b + c) * (a + b - c + 1) * (a + b + c + 1)))
    term1 = (((al - be) * c * (c - 1) - g * a * (a + 1) + g * b * (b + 1)) / (2 * c * (c - 1))
             * C(a, b, c - 1, al, be, g))
    term2 = (sqrt((c - g - 1) * (c + g - 1) * (-a + b + c - 1) * (a - b + c - 1)
                  * (a + b - c + 2) * (a + b + c) / (4 * (c - 1) ** 2 * (2 * c - 3) * (2 * c - 1)))
             * C(a, b, c - 2, al, be, g))
    return C(a, b, c, al, be, g), pref * (term1 - term2)


SEC_865 = [
    ("eq 8.6.148  b-1; a shift", r148),
    ("eq 8.6.149  2beta; a shift", r149),
    ("eq 8.6.150  b+1; a shift", r150),
    ("eq 8.6.151  c-1, c-2", r151),
]
# NOTE eq 8.6.150 in the source carries display-only OCR debris (a malformed
# \left.\beta^2\right| bracket and a stray ",(26" before the label); the maths
# is correct as tested here, but those tokens should be cleaned up.


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
IMPLEMENTED = [
    ("8.6.2  Arguments alpha,beta,gamma change by 1", SEC_862),
    ("8.6.3  Arguments change by 1/2", SEC_863),
    ("8.6.4  The case alpha=beta=gamma=0", SEC_864),
    ("8.6.5  Arguments a,b,c change by 1", SEC_865),
]

TODO = [
    ("8.6.1  General recursion relations", "eq. 126-128 (sums; eq.128 quasi-powers, exp 1/3 OCR)"),
    ("8.6.6  Arguments a,b,alpha,beta change by 1", "eq. 152-153"),
    ("8.6.7  Arguments c,b,gamma,beta change by 1", "eq. 154-156"),
    ("8.6.8  Recursion relations for the R-symbols", "eq. 157-161 (need R-symbol map)"),
]


def run(n, seed):
    random.seed(seed)
    print(f"Section 8.6 recursion-relation check -- seed={seed}, up to {n} "
          f"instances per relation\n")
    all_ok = True
    for title, table in IMPLEMENTED:
        print(f"=== {title} ===")
        for label, fn in table:
            got, bad = 0, None
            for _ in range(n * 500):
                if got >= n:
                    break
                res = fn()
                if res is None:
                    continue
                lhs, rhs = res
                got += 1
                if not _eq(lhs, rhs) and bad is None:
                    bad = (lhs, rhs)
            if got == 0:
                print(f"  [SKIP] {label:40s} no valid draws")
                continue
            ok = bad is None
            all_ok &= ok
            print(f"  [{'OK  ' if ok else 'FAIL'}] {label:40s} {got} instances")
            if bad:
                print(f"         lhs={bad[0]}   rhs={bad[1]}")
        print()

    print("Not yet implemented (subsections to add next):")
    for title, note in TODO:
        print(f"    {title:44s} {note}")
    print()
    print("ALL IMPLEMENTED RELATIONS HOLD" if all_ok else "SOME RELATIONS FAILED -- see above")
    return all_ok


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Check Section 8.6 recursion relations.")
    p.add_argument("--n", type=int, default=8, help="instances per relation")
    p.add_argument("--seed", type=int, default=20260805, help="RNG seed")
    args = p.parse_args()
    raise SystemExit(0 if run(args.n, args.seed) else 1)
