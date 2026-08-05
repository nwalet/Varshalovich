#!/usr/bin/env python3
r"""
Numerical check of the Clebsch-Gordan recursion relations listed in Section 8.6
of Varshalovich, Moskalev & Khersonskii, "Quantum Theory of Angular Momentum".

Each relation is an identity among Clebsch-Gordan coefficients whose arguments
are shifted by 1, 1/2, ... .  We draw random valid angular momenta/projections,
evaluate both sides with sympy's clebsch_gordan (coefficients that fall outside
the physical domain are treated as 0), and compare.

The tests are SPLIT BY SUBSECTION, matching the book:
    8.6.1  General recursion relations          eq. 126-128   [done]
    8.6.2  Arguments alpha,beta,gamma change 1  eq. 129-133   [done]
    8.6.3  Arguments change by 1/2              eq. 134-143   [done]
    8.6.4  The case alpha=beta=gamma=0          eq. 144-147   [done]
    8.6.5  Arguments a,b,c change by 1          eq. 148-151   [done]
    8.6.6  Arguments a,b,alpha,beta change 1    eq. 152-153   [done]
    8.6.7  Arguments c,b,gamma,beta change 1    eq. 154-156   [done]
    8.6.8  Recursion relations for R-symbols    eq. 157-161   [done, 158 skipped]

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
from sympy.physics.wigner import clebsch_gordan, wigner_3j

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


def qp(x, n):
    """quasi-power  x^(n) = x!/(x-n)!  (used in eq. 128)."""
    return fac(x) / fac(x - n)


# ===========================================================================
# 8.6.1  General recursion relations                 (eq. 126-128)
# ===========================================================================
def r126():
    # Yutsis-Bandzaitis: sum over c' of C_{a al, b-k, be-k}^{c', ga-k},
    # parametrised by k (integer or half-integer, 0 < k <= (b+be)/2).
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    g = al + be
    kmax = (b + be) / 2
    ks = [Rational(m, 2) for m in range(1, int(2 * kmax) + 1) if Rational(m, 2) <= kmax]
    if not ks:
        return None
    k = random.choice(ks)
    pref = sqrt(fac(b + be - 2 * k) * fac(c + g) * fac(a + b - c) * fac(-a + b + c)
                * fac(a + b + c + 1) * (2 * c + 1) / (fac(b + be) * fac(c - g) * fac(a - b + c)))
    tot = S.Zero
    for i in range(int(2 * k) + 1):
        cp = c - k + i
        cv = C(a, b - k, cp, al, be - k, g - k)
        if cv == 0:
            continue
        tot += ((-1) ** int(cp + k - c) * cv * fac(cp - k + c) * fac(2 * k)
                / (fac(c + k - cp) * fac(c + cp + k + 1) * fac(cp + k - c))
                * sqrt(fac(cp + k - g) * fac(a - b + cp + k) * (2 * cp + 1)
                       / (fac(cp - k + g) * fac(-a + b + cp - k)
                          * fac(a + b - cp - k) * fac(a + b + cp - k + 1))))
    return C(a, b, c, al, be, g), pref * tot


def r127():
    # Stone: sum over b' of C_{a-b, al-be, b', 2be}^{c ga}
    # (needs a-b >= |al-be| >= 0; b' integer, b'+2b even, |2be|,|-a+b+c| <= b' <= 2b, a-b+c)
    b = rj(HALF, 2); a = rj(b, 3); c = rc(a, b)
    al = rproj(a); be = rproj(b); g = al + be
    if abs(g) > c or a - b < abs(al - be) or C(a, b, c, al, be, g) == 0:
        return None
    pref = sqrt(fac(a + b + c + 1) * fac(a - b + c) * fac(a + b - c) * fac(b + be) * fac(b - be)
                * fac(a - b + al - be) * fac(a - b - al + be)
                / (fac(-a + b + c) * fac(a + al) * fac(a - al)))
    lo = max(abs(2 * be), abs(-a + b + c)); hi = min(2 * b, a - b + c)
    tot = S.Zero; nterm = 0
    bp = lo
    while bp <= hi:
        if Integer(bp + 2 * b) % 2 != 0:
            bp += 1
            continue
        cv = C(a - b, bp, c, al - be, 2 * be, g)
        if cv != 0:
            tot += ((-1) ** int(b - bp / 2) * cv * (2 * bp + 1) * fac(b + bp / 2)
                    / (fac(2 * b + bp + 1) * fac(b - bp / 2) * fac(bp / 2 + be) * fac(bp / 2 - be))
                    * sqrt(fac(-a + b + c + bp) * fac(bp + 2 * be) * fac(bp - 2 * be)
                           / (fac(a - b + c + bp + 1) * fac(a - b + c - bp) * fac(a - b - c + bp))))
            nterm += 1
        bp += 1
    if nterm == 0:
        return None
    return C(a, b, c, al, be, g), pref * tot


def r128():
    # quasi-power form of eq.126 (well-defined only for b+be >= 2k; the source
    # k-range 1/2 <= k <= (b-kappa)/2 must be intersected with that).
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    g = al + be
    kappa = 0 if b.is_integer else HALF
    kmax = (b - kappa) / 2
    ks = [Rational(m, 2) for m in range(1, int(2 * kmax) + 1)
          if Rational(m, 2) <= kmax and b + be >= m]
    if not ks:
        return None
    k = random.choice(ks)
    tot = S.Zero
    for j in range(int(2 * k) + 1):
        kp = -k + j
        cv = C(a, b - k, c + kp, al, be - k, g - k)
        if cv == 0:
            continue
        inner = (qp(c + g, k - kp) * qp(c - g + k + kp, k + kp) * qp(a + b - c, k + kp)
                 * qp(-a + b + c, k - kp) * qp(a + b + c + 1, k - kp)
                 * qp(a - b + c + k + kp, k + kp) * (2 * c + 2 * kp + 1))
        tot += ((-1) ** int(k + kp) * cv * sqrt(inner) * qp(2 * k, k + kp)
                / (qp(2 * c + k + kp + 1, 2 * k + 1) * qp(k + kp, k + kp)))
    return C(a, b, c, al, be, g), sqrt((2 * c + 1) / qp(b + be, 2 * k)) * tot


SEC_861 = [
    ("eq 8.6.126  Yutsis-Bandzaitis (sum c')", r126),
    ("eq 8.6.127  Stone (sum b')", r127),
    ("eq 8.6.128  quasi-power form", r128),
]


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


# ===========================================================================
# 8.6.6  Arguments a, b, alpha, beta change by 1      (eq. 152-153 + a middle
#        gather).  a shifted by -1,0,+1; alpha-/+1, beta+/-1 (so alpha+beta
#        stays gamma).
# ===========================================================================
def r152():
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    if b < 1:
        return None
    s = random.choice([1, -1]); g = al + be
    lhs = sqrt((b + s * be) * (b + s * be + 1)) * C(a, b - 1, c, al, be, g)
    rhs = (1 / (2 * a * (2 * a + 1))
           * sqrt((a + s * al) * (a + s * al - 1) * (-a + b + c) * (-a + b + c + 1) * (a - b + c) * (a - b + c + 1))
           * C(a - 1, b, c, al - s, be + s, g)
           - s / (2 * a * (a + 1))
           * sqrt((a + s * al) * (a - s * al + 1) * (-a + b + c) * (a - b + c + 1) * (a + b - c) * (a + b + c + 1))
           * C(a, b, c, al - s, be + s, g)
           + 1 / (2 * (a + 1) * (2 * a + 1))
           * sqrt((a - s * al + 1) * (a - s * al + 2) * (a + b - c) * (a + b - c + 1) * (a + b + c + 1) * (a + b + c + 2))
           * C(a + 1, b, c, al - s, be + s, g))
    return lhs, rhs


def r152b():
    # the middle gather (labelled chap8:eq:152b): 2beta-type, no a/b shift.
    # The {a(a+1)+b(b+1)-c(c+1)} factor sits OUTSIDE the radical.
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    s = random.choice([1, -1]); g = al + be
    lhs = sqrt((b - s * be) * (b + s * be + 1)) * C(a, b, c, al, be, g)
    rhs = (s / (2 * a * (2 * a + 1))
           * sqrt((a + s * al) * (a + s * al - 1) * (-a + b + c + 1) * (a - b + c) * (a + b - c) * (a + b + c + 1))
           * C(a - 1, b, c, al - s, be + s, g)
           - 1 / (2 * a * (a + 1)) * (a * (a + 1) + b * (b + 1) - c * (c + 1))
           * sqrt((a + s * al) * (a - s * al + 1))
           * C(a, b, c, al - s, be + s, g)
           - s / (2 * (a + 1) * (2 * a + 1))
           * sqrt((a - s * al + 1) * (a - s * al + 2) * (-a + b + c) * (a - b + c + 1) * (a + b - c + 1) * (a + b + c + 2))
           * C(a + 1, b, c, al - s, be + s, g))
    return lhs, rhs


def r153():
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    s = random.choice([1, -1]); g = al + be
    lhs = sqrt((b - s * be) * (b - s * be + 1)) * C(a, b + 1, c, al, be, g)
    rhs = (1 / (2 * a * (2 * a + 1))
           * sqrt((a + s * al) * (a + s * al - 1) * (a + b - c) * (a + b - c + 1) * (a + b + c + 1) * (a + b + c + 2))
           * C(a - 1, b, c, al - s, be + s, g)
           + s / (2 * a * (a + 1))
           * sqrt((a + s * al) * (a - s * al + 1) * (-a + b + c + 1) * (a - b + c) * (a + b - c + 1) * (a + b + c + 2))
           * C(a, b, c, al - s, be + s, g)
           + 1 / (2 * (a + 1) * (2 * a + 1))
           * sqrt((a - s * al + 1) * (a - s * al + 2) * (-a + b + c) * (-a + b + c + 1) * (a - b + c) * (a - b + c + 1))
           * C(a + 1, b, c, al - s, be + s, g))
    return lhs, rhs


SEC_866 = [
    ("eq 8.6.152  b-1; a shift, al-/+1", r152),
    ("eq 8.6.152b middle gather", r152b),
    ("eq 8.6.153  b+1; a shift, al-/+1", r153),
]


# ===========================================================================
# 8.6.7  Arguments c, b, gamma, beta change by 1      (eq. 154-156 + a middle
#        gather).  c shifted by -1,0,+1; beta+/-1, gamma+/-1.
# ===========================================================================
def r154():
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    if c < 1:
        return None
    g = al + be
    lhs = sqrt((-a + b + c) * (a - b + c) * (a + b - c + 1) * (a + b + c + 1) * (2 * c - 1) / (2 * c + 1)) * C(a, b, c, al, be, g)
    rhs = (sqrt((b + be) * (b - be + 1) * (c + g) * (c + g - 1)) * C(a, b, c - 1, al, be - 1, g - 1)
           - 2 * be * sqrt(c ** 2 - g ** 2) * C(a, b, c - 1, al, be, g)
           - sqrt((b - be) * (b + be + 1) * (c - g) * (c - g - 1)) * C(a, b, c - 1, al, be + 1, g + 1))
    return lhs, rhs


def r155():
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    if b < 1 or c < 1:
        return None
    s = random.choice([1, -1]); g = al + be
    lhs = sqrt((b - s * be) * (b - s * be + 1)) * C(a, b - 1, c, al, be, g)
    rhs = (1 / (2 * c)
           * sqrt((c + s * g) * (c + s * g - 1) * (a - b + c) * (a - b + c + 1) * (a + b - c) * (a + b - c + 1) / ((2 * c - 1) * (2 * c + 1)))
           * C(a, b, c - 1, al, be - s, g - s)
           + s / (2 * c * (c + 1))
           * sqrt((c + s * g) * (c - s * g + 1) * (-a + b + c) * (a - b + c + 1) * (a + b - c) * (a + b + c + 1))
           * C(a, b, c, al, be - s, g - s)
           + 1 / (2 * (c + 1))
           * sqrt((c - s * g + 1) * (c - s * g + 2) * (-a + b + c) * (-a + b + c + 1) * (a + b + c + 1) * (a + b + c + 2) / ((2 * c + 1) * (2 * c + 3)))
           * C(a, b, c + 1, al, be - s, g - s))
    return lhs, rhs


def r155b():
    # the middle gather (labelled chap8:eq:155b): 2beta-type, no b/c shift.
    # The {-a(a+1)+b(b+1)+c(c+1)} factor sits OUTSIDE the radical.
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    if c < 1:
        return None
    s = random.choice([1, -1]); g = al + be
    lhs = sqrt((b + s * be) * (b - s * be + 1)) * C(a, b, c, al, be, g)
    rhs = (s / (2 * c)
           * sqrt((c + s * g) * (c + s * g - 1) * (-a + b + c) * (a - b + c) * (a + b - c + 1) * (a + b + c + 1) / ((2 * c - 1) * (2 * c + 1)))
           * C(a, b, c - 1, al, be - s, g - s)
           + 1 / (2 * c * (c + 1)) * (-a * (a + 1) + b * (b + 1) + c * (c + 1))
           * sqrt((c + s * g) * (c - s * g + 1))
           * C(a, b, c, al, be - s, g - s)
           - s / (2 * (c + 1))
           * sqrt((c - s * g + 1) * (c - s * g + 2) * (-a + b + c + 1) * (a - b + c + 1) * (a + b - c) * (a + b + c + 2) / ((2 * c + 1) * (2 * c + 3)))
           * C(a, b, c + 1, al, be - s, g - s))
    return lhs, rhs


def r156():
    cfg = cfg_master()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    if c < 1:
        return None
    s = random.choice([1, -1]); g = al + be
    lhs = sqrt((b + s * be) * (b + s * be + 1)) * C(a, b + 1, c, al, be, g)
    rhs = (1 / (2 * c)
           * sqrt((c + s * g) * (c + s * g - 1) * (-a + b + c) * (-a + b + c + 1) * (a + b + c + 1) * (a + b + c + 2) / ((2 * c - 1) * (2 * c + 1)))
           * C(a, b, c - 1, al, be - s, g - s)
           - s / (2 * c * (c + 1))
           * sqrt((c + s * g) * (c - s * g + 1) * (-a + b + c + 1) * (a - b + c) * (a + b - c + 1) * (a + b + c + 2))
           * C(a, b, c, al, be - s, g - s)
           + 1 / (2 * (c + 1))
           * sqrt((c - s * g + 1) * (c - s * g + 2) * (a - b + c) * (a - b + c + 1) * (a + b - c) * (a + b - c + 1) / ((2 * c + 1) * (2 * c + 3)))
           * C(a, b, c + 1, al, be - s, g - s))
    return lhs, rhs


SEC_867 = [
    ("eq 8.6.154  c-1; beta,gamma shift", r154),
    ("eq 8.6.155  b-1; c shift", r155),
    ("eq 8.6.155b middle gather", r155b),
    ("eq 8.6.156  b+1; c shift", r156),
]


# ===========================================================================
# 8.6.8  Recursion relations for the Regge R-symbols  (eq. 157-161)
#
# The Regge R-symbol ||R|| equals the 3jm symbol (Sec. 8.1.3, eq. 13-15):
#   R = [[-a+b+c, a-b+c, a+b-c],
#        [ a+al ,  b+be,  c+ga ],
#        [ a-al ,  b-be,  c-ga ]]     (alpha+beta+gamma = 0).
# Each relation shifts individual entries by +/-1 (line sums move J -> J+/-1)
# and the shifted array is evaluated as its own 3jm symbol.
# ===========================================================================
def Rmat(a, b, c, al, be, ga):
    return [[-a + b + c, a - b + c, a + b - c],
            [a + al, b + be, c + ga],
            [a - al, b - be, c - ga]]


def shift(R, *changes):
    m = [row[:] for row in R]
    for (i, j, d) in changes:
        m[i][j] = m[i][j] + d
    return m


def Rval(m):
    """value of the Regge symbol ||m|| = 3jm recovered from the array
    (0 unless it is a valid, line-balanced non-negative-integer array)."""
    for row in m:
        for e in row:
            if e < 0:
                return S.Zero
    a = Rational(m[1][0] + m[2][0], 2); b = Rational(m[1][1] + m[2][1], 2); c = Rational(m[1][2] + m[2][2], 2)
    al = Rational(m[1][0] - m[2][0], 2); be = Rational(m[1][1] - m[2][1], 2); ga = Rational(m[1][2] - m[2][2], 2)
    if al + be + ga != 0 or c < abs(a - b) or c > a + b or not (a + b + c).is_integer:
        return S.Zero
    if abs(al) > a or abs(be) > b or abs(ga) > c:
        return S.Zero
    return wigner_3j(a, b, c, al, be, ga)


def cfg_3jm(jmax=3):
    """valid (a,b,c,alpha,beta) with gamma=-(alpha+beta) and a non-zero 3jm."""
    a = rj(HALF, jmax); b = rj(HALF, jmax); c = rc(a, b)
    al = rproj(a); be = rproj(b); ga = -(al + be)
    if abs(ga) > c or Rval(Rmat(a, b, c, al, be, ga)) == 0:
        return None
    return (a, b, c, al, be)


def _rsym(build):
    """shared driver: build the list of (signed) terms, require the relation to
    be non-vacuous (>= 2 terms non-zero), and return (sum, 0)."""
    cfg = cfg_3jm()
    if cfg is None:
        return None
    a, b, c, al, be = cfg
    R = Rmat(a, b, c, al, be, -(al + be))
    J = a + b + c
    terms = build(R, J)
    if sum(1 for t in terms if t != 0) < 2:
        return None
    return sum(terms, S.Zero), S.Zero


def r157():
    def build(R, J):
        return [sqrt(R[0][0] * (J + 1)) * Rval(R),
                -sqrt(R[1][1] * R[2][2]) * Rval(shift(R, (0, 0, -1), (1, 1, -1), (2, 2, -1))),
                sqrt(R[1][2] * R[2][1]) * Rval(shift(R, (0, 0, -1), (1, 2, -1), (2, 1, -1)))]
    return _rsym(build)


def r159():
    def build(R, J):
        return [sqrt(R[0][0] * R[1][1] * (J + 1)) * Rval(R),
                -(R[1][1] + R[1][2]) * sqrt(R[2][2]) * Rval(shift(R, (0, 0, -1), (1, 1, -1), (2, 2, -1))),
                -sqrt(R[1][2] * R[2][0] * (R[1][0] + 1))
                * Rval(shift(R, (0, 0, -1), (1, 0, 1), (1, 1, -1), (1, 2, -1), (2, 0, -1)))]
    return _rsym(build)


def r160():
    def build(R, J):
        return [(R[0][0] + R[0][1] + 1) * sqrt(R[1][1]) * Rval(R),
                sqrt(R[0][2] * (R[1][2] + 1) * (R[0][1] + 1))
                * Rval(shift(R, (0, 1, 1), (0, 2, -1), (1, 1, -1), (1, 2, 1))),
                -sqrt(R[0][0] * R[2][2] * (J + 1)) * Rval(shift(R, (0, 0, -1), (1, 1, -1), (2, 2, -1)))]
    return _rsym(build)


def r161():
    def build(R, J):
        return [(R[1][1] - R[2][2]) * sqrt(R[0][0] * R[1][2] * R[2][1])
                * Rval(shift(R, (0, 0, -1), (1, 2, -1), (2, 1, -1))),
                (R[2][2] - R[0][0]) * sqrt(R[0][2] * R[1][1] * R[2][0])
                * Rval(shift(R, (0, 2, -1), (1, 1, -1), (2, 0, -1))),
                (R[0][0] - R[1][1]) * sqrt(R[0][1] * R[1][0] * R[2][2])
                * Rval(shift(R, (0, 1, -1), (1, 0, -1), (2, 2, -1)))]
    return _rsym(build)


SEC_868 = [
    ("eq 8.6.157  Regge R-symbol", r157),
    ("eq 8.6.159  Regge R-symbol", r159),
    ("eq 8.6.160  Regge R-symbol", r160),
    ("eq 8.6.161  Regge R-symbol", r161),
]
# NOTE eq 8.6.158 is a "column" recursion: it shifts R_{1k}+1 and R_{2k}-1, so
# each array has equal COLUMN sums but row sums J+1, J-1, J.  A Regge symbol
# needs all six line sums equal, so every term vanishes as a standard 3jm and
# the relation is only trivially (vacuously) satisfied here.  Verifying it would
# require the generalised symbol of refs [86,103]; skipped.


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
IMPLEMENTED = [
    ("8.6.1  General recursion relations", SEC_861),
    ("8.6.2  Arguments alpha,beta,gamma change by 1", SEC_862),
    ("8.6.3  Arguments change by 1/2", SEC_863),
    ("8.6.4  The case alpha=beta=gamma=0", SEC_864),
    ("8.6.5  Arguments a,b,c change by 1", SEC_865),
    ("8.6.6  Arguments a,b,alpha,beta change by 1", SEC_866),
    ("8.6.7  Arguments c,b,gamma,beta change by 1", SEC_867),
    ("8.6.8  Recursion relations for the R-symbols", SEC_868),
]

TODO = [
    ("8.6.8  eq 158 only", "column recursion -- row-unbalanced arrays, vacuous as standard 3jm"),
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
