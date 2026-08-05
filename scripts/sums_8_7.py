#!/usr/bin/env python3
r"""
Numerical check of the sum rules for products of Clebsch-Gordan coefficients in
Section 8.7 of Varshalovich, Moskalev & Khersonskii, "Quantum Theory of Angular
Momentum".

Each identity fixes some external angular momenta/projections and sums a product
of Clebsch-Gordan coefficients over internal ones.  We pick random valid
externals, evaluate the sum with sympy, and compare with the closed-form RHS.
Tests are split by subsection.

Convention (book == sympy):  C_{a al, b be}^{c ga} == clebsch_gordan(a,b,c,al,be,ga)
Notation:  Pi(a,b,..) = [(2a+1)(2b+1)...]^{1/2}   (the book's \Fact{a b ...}).

Covered
    8.7.1  Sums involving one CG                 eq. 162-164   [done]
    8.7.2  Sums of products of two CG            eq. 165-172   [done]
    8.7.3  three CG (one 6j)                     eq. 173-180   [done]
    8.7.5  CG and one 6j                         eq. 192-199   [done]
    8.7.7  Additional sums of two CG             eq. 205-208   [done]

Not yet covered:
    8.7.4  four CG (one 9j)                      eq. 181-191   (some OCR damage)
    8.7.6  CG and one 9j                         eq. 200-204   (some OCR damage)

Usage:
    python3 sums_8_7.py [--n N] [--seed S]
"""
from __future__ import annotations

import argparse
import random

from sympy import Rational, Integer, sqrt, factorial as fac, S
from sympy.physics.wigner import clebsch_gordan, wigner_6j

HALF = Rational(1, 2)
TOL = S(10) ** (-18)


def _eq(x, y) -> bool:
    d = (x - y).evalf(30)
    return d.is_finite and d.is_real and abs(d) < TOL


def ph(n) -> int:
    return int(S.NegativeOne ** n)


def C(a, b, cc, al, be, ga):
    if a < 0 or b < 0 or cc < 0:
        return S.Zero
    if abs(al) > a or abs(be) > b or abs(ga) > cc:
        return S.Zero
    if cc < abs(a - b) or cc > a + b:
        return S.Zero
    if not (a - al).is_integer or not (b - be).is_integer or not (cc - ga).is_integer:
        return S.Zero
    if not (a + b + cc).is_integer or al + be != ga:
        return S.Zero
    return clebsch_gordan(a, b, cc, al, be, ga)


def Pi(*js):
    p = S.One
    for j in js:
        p *= (2 * j + 1)
    return sqrt(p)


def rj(lo, hi):
    return Rational(random.randint(int(2 * lo), int(2 * hi)), 2)


def rint(lo, hi):
    return Integer(random.randint(int(lo), int(hi)))


def proj(j):
    return [-j + i for i in range(int(2 * j) + 1)]


def crange(a, b):
    return [abs(a - b) + i for i in range(int(2 * min(a, b)) + 1)]


def w6(a, b, c, d, e, f):
    """Wigner 6j symbol {a b c; d e f} (0 outside the physical domain)."""
    try:
        return wigner_6j(a, b, c, d, e, f)
    except Exception:
        return S.Zero


def sgn(n):
    """(-1)^n kept symbolic (n may be a momentum sum, not always integer)."""
    return S.NegativeOne ** n


MR = [Rational(i, 2) for i in range(0, 9)]     # internal momentum range 0..4


# ===========================================================================
# 8.7.1  Sums involving one Clebsch-Gordan coefficient       (eq. 162-164)
# ===========================================================================
def r162():
    # sum_al C_{a al, b 0}^{a al} = Pi(a)^2 delta_{b0}
    a = rj(HALF, 3)
    b = rint(0, int(2 * a))                      # b integer, 0..2a (beta=0 needs b integer)
    lhs = sum((C(a, b, a, al, 0, al) for al in proj(a)), S.Zero)
    rhs = (2 * a + 1) * (1 if b == 0 else 0)
    return lhs, rhs


def r163():
    # sum_al (-1)^{a-al} C_{a al, a -al}^{c 0} = Pi(a) delta_{c0}
    a = rj(HALF, 3)
    c = rint(0, int(2 * a))
    lhs = sum((ph(a - al) * C(a, a, c, al, -al, 0) for al in proj(a)), S.Zero)
    rhs = sqrt(2 * a + 1) * (1 if c == 0 else 0)
    return lhs, rhs


def r164():
    # sum_{al be ga} (-1)^{c+ga} C_{a al,b be}^{c ga} [prod (j+-m)!]^{-1/2} = 0
    a = rj(HALF, 2); b = rj(HALF, 2); c = crange(a, b)[random.randrange(len(crange(a, b)))]
    tot = S.Zero
    for al in proj(a):
        for be in proj(b):
            ga = al + be
            cv = C(a, b, c, al, be, ga)
            if cv == 0:
                continue
            tot += (ph(c + ga) * cv
                    / sqrt(fac(a + al) * fac(a - al) * fac(b + be) * fac(b - be)
                           * fac(c + ga) * fac(c - ga)))
    return tot, S.Zero


SEC_871 = [
    ("eq 8.7.162  sum_al C_{a,b0}^{a}", r162),
    ("eq 8.7.163  sum_al (-1)^.. C_{a,a}^{c0}", r163),
    ("eq 8.7.164  weighted sum = 0", r164),
]


# ===========================================================================
# 8.7.2  Sums of products of two Clebsch-Gordan coefficients  (eq. 165-172)
# ===========================================================================
def r165():
    # sum_{al be} C_{a al,b be}^{c ga} C_{a al,b be}^{c' ga'} = d_{cc'} d_{ga ga'}
    a = rj(HALF, 3); b = rj(HALF, 3)
    cs = crange(a, b)
    c = random.choice(cs); cp = random.choice(cs)
    ga = random.choice(proj(c)); gap = random.choice(proj(cp))
    lhs = sum((C(a, b, c, al, be, ga) * C(a, b, cp, al, be, gap)
               for al in proj(a) for be in proj(b)), S.Zero)
    rhs = 1 if (c == cp and ga == gap) else 0
    return lhs, rhs


def r166():
    # sum_{al ga} C_{a al,b be}^{c ga} C_{a al,b' be'}^{c ga} = Pi(c)^2/Pi(b)^2 d_{bb'} d_{be be'}
    a = rj(HALF, 3); b = rj(HALF, 3)
    cs = crange(a, b)
    c = random.choice(cs)
    bp = random.choice([b] + crange(a, c))       # b' sometimes = b
    be = random.choice(proj(b)); bep = random.choice(proj(bp))
    lhs = sum((C(a, b, c, al, be, ga) * C(a, bp, c, al, bep, ga)
               for al in proj(a) for ga in proj(c)), S.Zero)
    rhs = (2 * c + 1) / (2 * b + 1) if (bp == b and be == bep) else 0
    return lhs, rhs


def r167():
    # sum_{al be} (-1)^{b+be} C_{a al,b be}^{c ga} C_{c' -ga', a al}^{b -be}
    #           = (-1)^{c+ga} Pi(b)/Pi(c) d_{cc'} d_{ga ga'}
    a = rj(HALF, 3); b = rj(HALF, 3)
    cs = crange(a, b)
    c = random.choice(cs); cp = random.choice(cs)
    ga = random.choice(proj(c)); gap = random.choice(proj(cp))
    lhs = sum((ph(b + be) * C(a, b, c, al, be, ga) * C(cp, a, b, -gap, al, -be)
               for al in proj(a) for be in proj(b)), S.Zero)
    rhs = ph(c + ga) * Pi(b) / Pi(c) if (c == cp and ga == gap) else 0
    return lhs, rhs


def r168():
    # sum_{al be} (-1)^{a+al} C_{b be,a al}^{c ga} C_{a -al,c' ga'}^{b be}
    #           = Pi(b)/Pi(c) d_{cc'} d_{ga ga'}
    a = rj(HALF, 3); b = rj(HALF, 3)
    cs = crange(a, b)
    c = random.choice(cs); cp = random.choice(cs)
    ga = random.choice(proj(c)); gap = random.choice(proj(cp))
    lhs = sum((ph(a + al) * C(b, a, c, be, al, ga) * C(a, cp, b, -al, gap, be)
               for al in proj(a) for be in proj(b)), S.Zero)
    rhs = Pi(b) / Pi(c) if (c == cp and ga == gap) else 0
    return lhs, rhs


def r169():
    # sum_{al be} C_{c ga,b be}^{a al} C_{c' -ga',a al}^{b be}
    #           = (-1)^{b-a-ga} Pi(a b)/Pi(c)^2 d_{cc'} d_{ga ga'}
    a = rj(HALF, 3); b = rj(HALF, 3)
    cs = crange(a, b)
    c = random.choice(cs); cp = random.choice(cs)
    ga = random.choice(proj(c)); gap = random.choice(proj(cp))
    lhs = sum((C(c, b, a, ga, be, al) * C(cp, a, b, -gap, al, be)
               for al in proj(a) for be in proj(b)), S.Zero)
    rhs = ph(b - a - ga) * Pi(a, b) / (2 * c + 1) if (c == cp and ga == gap) else 0
    return lhs, rhs


def r170():
    # sum_{c ga} C_{a al,b be}^{c ga} C_{a al',b be'}^{c ga} = d_{al al'} d_{be be'}
    a = rj(HALF, 3); b = rj(HALF, 3)
    al = random.choice(proj(a)); alp = random.choice(proj(a))
    be = random.choice(proj(b)); bep = random.choice(proj(b))
    lhs = sum((C(a, b, c, al, be, al + be) * C(a, b, c, alp, bep, alp + bep)
               for c in crange(a, b) for ga in proj(c)
               if ga == al + be and ga == alp + bep), S.Zero)
    rhs = 1 if (al == alp and be == bep) else 0
    return lhs, rhs


def r171():
    # sum_{c ga} Pi(c)^2 C_{a al,c ga}^{b be} C_{a al',c ga}^{b be'} = Pi(b)^2 d_{al al'} d_{be be'}
    a = rj(HALF, 3); b = rj(HALF, 3)
    al = random.choice(proj(a)); alp = random.choice(proj(a))
    be = random.choice(proj(b)); bep = random.choice(proj(b))
    lhs = sum(((2 * c + 1) * C(a, c, b, al, ga, be) * C(a, c, b, alp, ga, bep)
               for c in crange(a, b) for ga in proj(c)), S.Zero)
    rhs = (2 * b + 1) if (al == alp and be == bep) else 0
    return lhs, rhs


def r172():
    # sum_{a al} (-1)^{a-al} Pi(a)^2 C_{a al,b be}^{c ga} C_{a al,c ga'}^{b be'}
    #          = Pi(b c) d_{ga,-ga'} d_{be,-be'}
    b = rj(HALF, 2); c = rj(HALF, 2)
    be = random.choice(proj(b)); bep = random.choice(proj(b))
    ga = random.choice(proj(c)); gap = random.choice(proj(c))
    al = ga - be                                  # first CG pins alpha = gamma - beta
    lhs = sum(((-1) ** int(a - al) * (2 * a + 1)
               * C(a, b, c, al, be, ga) * C(a, c, b, al, gap, bep)
               for a in crange(b, c)), S.Zero)
    rhs = Pi(b, c) if (ga == -gap and be == -bep) else 0
    return lhs, rhs


SEC_872 = [
    ("eq 8.7.165  orthogonality (c,ga)", r165),
    ("eq 8.7.166  orthogonality (b,be)", r166),
    ("eq 8.7.167  recoupled orthogonality", r167),
    ("eq 8.7.168  recoupled orthogonality", r168),
    ("eq 8.7.169  recoupled orthogonality", r169),
    ("eq 8.7.170  completeness (al,be)", r170),
    ("eq 8.7.171  weighted completeness", r171),
    ("eq 8.7.172  weighted sum over a", r172),
]


# ===========================================================================
# 8.7.3  Sums of products of three CG (one 6j)               (eq. 173-180)
#
# All have the form   sum_{al be de} [phase] C C C
#                        = kappa Pi(..) C_{c ga, f phi}^{e eps} {a b c; e f d}
# with kappa1 = (-1)^{b+c+d+f}, kappa2 = (-1)^{a+b+e+f}, and eps = ga + phi.
# ===========================================================================
def _cfg873():
    """(a,b,c,d,e,f, ga, phi, eps) with {a b c; e f d} != 0 and C_{c,f}^{e} != 0."""
    a = b = c = d = e = f = None
    for _ in range(400):
        a, b, c, d, e, f = (rj(HALF, 2) for _ in range(6))
        if w6(a, b, c, e, f, d) != 0:
            break
    else:
        return None
    ga = random.choice(proj(c)); phi = random.choice(proj(f)); eps = ga + phi
    if abs(eps) > e or C(c, f, e, ga, phi, eps) == 0:
        return None
    return (a, b, c, d, e, f, ga, phi, eps)


def _triple(a, b, d, term):
    return sum((term(al, be, de) for al in proj(a) for be in proj(b) for de in proj(d)), S.Zero)


def r173():
    P = _cfg873()
    if P is None:
        return None
    a, b, c, d, e, f, ga, phi, eps = P
    lhs = _triple(a, b, d, lambda al, be, de:
                  C(a, b, c, al, be, ga) * C(d, b, e, de, be, eps) * C(a, f, d, al, phi, de))
    rhs = ph(b + c + d + f) * Pi(c, d) * C(c, f, e, ga, phi, eps) * w6(a, b, c, e, f, d)
    return lhs, rhs


def r174():
    P = _cfg873()
    if P is None:
        return None
    a, b, c, d, e, f, ga, phi, eps = P
    lhs = _triple(a, b, d, lambda al, be, de:
                  C(b, c, a, be, ga, al) * C(b, e, d, be, eps, de) * C(a, f, d, al, phi, de))
    rhs = ph(b + c + d + f) * Pi(a, d, d) / Pi(e) * C(c, f, e, ga, phi, eps) * w6(a, b, c, e, f, d)
    return lhs, rhs


def r175():
    P = _cfg873()
    if P is None:
        return None
    a, b, c, d, e, f, ga, phi, eps = P
    lhs = _triple(a, b, d, lambda al, be, de:
                  C(b, a, c, be, al, ga) * C(b, d, e, be, de, eps) * C(a, f, d, al, phi, de))
    rhs = ph(a + b + e + f) * Pi(c, d) * C(c, f, e, ga, phi, eps) * w6(a, b, c, e, f, d)
    return lhs, rhs


def r176():
    P = _cfg873()
    if P is None:
        return None
    a, b, c, d, e, f, ga, phi, eps = P
    lhs = _triple(a, b, d, lambda al, be, de: ph(a - al) *
                  C(a, b, c, al, be, ga) * C(d, b, e, de, be, eps) * C(d, a, f, de, -al, phi))
    rhs = ph(b + c + d + f) * Pi(c, f) * C(c, f, e, ga, phi, eps) * w6(a, b, c, e, f, d)
    return lhs, rhs


def r177():
    P = _cfg873()
    if P is None:
        return None
    a, b, c, d, e, f, ga, phi, eps = P
    lhs = _triple(a, b, d, lambda al, be, de: ph(b + be) *
                  C(a, b, c, al, be, ga) * C(b, e, d, -be, eps, de) * C(a, f, d, al, phi, de))
    rhs = ph(b + c + d + f) * Pi(c, d, d) / Pi(e) * C(c, f, e, ga, phi, eps) * w6(a, b, c, e, f, d)
    return lhs, rhs


def r178():
    P = _cfg873()
    if P is None:
        return None
    a, b, c, d, e, f, ga, phi, eps = P
    lhs = _triple(a, b, d, lambda al, be, de: ph(a - al) *
                  C(b, a, c, be, al, ga) * C(b, d, e, be, de, eps) * C(d, a, f, de, -al, phi))
    rhs = ph(a + b + e + f) * Pi(c, f) * C(c, f, e, ga, phi, eps) * w6(a, b, c, e, f, d)
    return lhs, rhs


def r179():
    P = _cfg873()
    if P is None:
        return None
    a, b, c, d, e, f, ga, phi, eps = P
    lhs = _triple(a, b, d, lambda al, be, de: ph(b + be) *
                  C(b, c, a, -be, ga, al) * C(d, b, e, de, be, eps) * C(a, f, d, al, phi, de))
    rhs = ph(b + c + d + f) * Pi(a, d) * C(c, f, e, ga, phi, eps) * w6(a, b, c, e, f, d)
    return lhs, rhs


def r180():
    P = _cfg873()
    if P is None:
        return None
    a, b, c, d, e, f, ga, phi, eps = P
    lhs = _triple(a, b, d, lambda al, be, de: ph(a - al) *
                  C(a, c, b, al, -ga, be) * C(d, b, e, de, -be, eps) * C(a, f, d, al, phi, de))
    rhs = ph(b + c + d + f) * Pi(b, d) * C(c, f, e, ga, phi, eps) * w6(a, b, c, e, f, d)
    return lhs, rhs


SEC_873 = [
    ("eq 8.7.173  three CG + 6j", r173),
    ("eq 8.7.174  three CG + 6j", r174),
    ("eq 8.7.175  three CG + 6j", r175),
    ("eq 8.7.176  three CG + 6j", r176),
    ("eq 8.7.177  three CG + 6j", r177),
    ("eq 8.7.178  three CG + 6j", r178),
    ("eq 8.7.179  three CG + 6j", r179),
    ("eq 8.7.180  three CG + 6j", r180),
]


# ===========================================================================
# 8.7.5  Sums of products of CG and one 6j symbol            (eq. 192-199)
#
# Form  sum_{internal} [phase] Pi(..) C C {6j} = C C   (two CG on the RHS).
# The internal (summed) momentum differs per identity; its projection is
# pinned by the coefficients, so only the momentum is summed (over MR).
# ===========================================================================
def r192():
    a, b, c, d, f = (rj(HALF, 2) for _ in range(5))
    al = random.choice(proj(a)); be = random.choice(proj(b)); phi = random.choice(proj(f))
    ga = al + be; de = al + phi; eps = al + be + phi
    if abs(ga) > c or abs(de) > d:
        return None
    rhs = C(a, b, c, al, be, ga) * C(a, f, d, al, phi, de)
    if rhs == 0:
        return None
    lhs = sum((sgn(2 * e) * Pi(c, d) * C(b, d, e, be, de, eps) * C(f, c, e, phi, ga, eps)
               * w6(a, b, c, e, f, d) for e in MR), S.Zero)
    return lhs, rhs


def r193():
    a, b, c, d, e = (rj(HALF, 2) for _ in range(5))
    al = random.choice(proj(a)); be = random.choice(proj(b)); de = random.choice(proj(d))
    ga = al + be; eps = de + be; phi = de + al + be
    if abs(ga) > c or abs(eps) > e:
        return None
    rhs = C(a, b, c, al, be, ga) * C(d, b, e, de, be, eps)
    if rhs == 0:
        return None
    lhs = sum((sgn(c + d + f) * Pi(c, e) * C(e, a, f, eps, al, phi) * C(d, c, f, de, ga, phi)
               * w6(b, a, c, f, d, e) for f in MR), S.Zero)
    return lhs, rhs


def r194():
    a, b, d, e, f = (rj(HALF, 2) for _ in range(5))
    be = random.choice(proj(b)); de = random.choice(proj(d)); phi = random.choice(proj(f))
    eps = be + de; al = phi + de; ga = be - phi
    if abs(eps) > e or abs(al) > a:
        return None
    rhs = C(b, d, e, be, de, eps) * C(f, d, a, phi, de, al)
    if rhs == 0:
        return None
    lhs = sum((sgn(2 * e - d + al + phi) * Pi(a, e) * C(f, b, c, -phi, be, ga)
               * C(e, a, c, eps, -al, ga) * w6(c, f, b, d, e, a) for c in MR), S.Zero)
    return lhs, rhs


def r195():
    a, b, d, e, f = (rj(HALF, 2) for _ in range(5))
    al = random.choice(proj(a)); be = random.choice(proj(b)); phi = random.choice(proj(f))
    de = al + phi; eps = al + be + phi; ga = al + be
    if abs(de) > d or abs(eps) > e:
        return None
    rhs = C(a, f, d, al, phi, de) * C(b, e, d, -be, eps, de)
    if rhs == 0:
        return None
    lhs = sum((sgn(c + d - be - phi) * (2 * d + 1) * C(a, b, c, al, be, ga)
               * C(f, e, c, -phi, eps, ga) * w6(a, b, c, e, f, d) for c in MR), S.Zero)
    return lhs, rhs


def r196():
    a, b, d, e, f = (rj(HALF, 2) for _ in range(5))
    al = random.choice(proj(a)); be = random.choice(proj(b)); phi = random.choice(proj(f))
    de = al + phi; eps = al + be + phi; ga = al + be
    if abs(de) > d or abs(eps) > e:
        return None
    rhs = C(b, d, e, be, de, eps) * C(a, f, d, al, phi, de)
    if rhs == 0:
        return None
    lhs = sum((sgn(2 * e) * Pi(c, d) * C(a, b, c, al, be, ga) * C(f, c, e, phi, ga, eps)
               * w6(a, b, c, e, f, d) for c in MR), S.Zero)
    return lhs, rhs


def r197():
    a, b, c, d, e = (rj(HALF, 2) for _ in range(5))
    al = random.choice(proj(a)); be = random.choice(proj(b)); de = random.choice(proj(d))
    eps = be + al; ga = de + eps; phi = be + de
    if abs(eps) > e or abs(ga) > c:
        return None
    rhs = C(b, a, e, be, al, eps) * C(d, e, c, de, eps, ga)
    if rhs == 0:
        return None
    lhs = sum((sgn(2 * c) * Pi(e, f) * C(b, d, f, be, de, phi) * C(a, f, c, al, phi, ga)
               * w6(a, b, e, d, c, f) for f in MR), S.Zero)
    return lhs, rhs


def r198():
    a, b, d, e, f = (rj(HALF, 2) for _ in range(5))
    al = random.choice(proj(a)); be = random.choice(proj(b)); phi = random.choice(proj(f))
    de = al + phi; eps = al + be + phi; ga = al + be
    if abs(de) > d or abs(eps) > e:
        return None
    rhs = C(b, e, d, -be, eps, de) * C(a, f, d, al, phi, de)
    if rhs == 0:
        return None
    lhs = sum((sgn(e + d - be) * sqrt((2 * c + 1) * (2 * d + 1) ** 2 / (2 * e + 1))
               * C(a, b, c, al, be, ga) * C(f, c, e, phi, ga, eps) * w6(a, b, c, e, f, d)
               for c in MR), S.Zero)
    return lhs, rhs


def r199():
    a, b, d, e, f = (rj(HALF, 2) for _ in range(5))
    be = random.choice(proj(b)); de = random.choice(proj(d)); phi = random.choice(proj(f))
    eps = de + be; al = phi + de; ga = be - phi
    if abs(eps) > e or abs(al) > a:
        return None
    rhs = C(d, b, e, de, be, eps) * C(f, d, a, phi, de, al)
    if rhs == 0:
        return None
    lhs = sum((sgn(2 * e) * sqrt((2 * a + 1) * (2 * c + 1) ** 2 / (2 * b + 1))
               * C(f, c, b, phi, ga, be) * C(c, a, e, ga, al, eps) * w6(c, f, b, d, e, a)
               for c in MR), S.Zero)
    return lhs, rhs


SEC_875 = [
    ("eq 8.7.192  CG + 6j (sum e)", r192),
    ("eq 8.7.193  CG + 6j (sum f)", r193),
    ("eq 8.7.194  CG + 6j (sum c)", r194),
    ("eq 8.7.195  CG + 6j (sum c)", r195),
    ("eq 8.7.196  CG + 6j (sum c)", r196),
    ("eq 8.7.197  CG + 6j (sum f)", r197),
    ("eq 8.7.198  CG + 6j (sum c)", r198),
    ("eq 8.7.199  CG + 6j (sum c)", r199),
]


# ===========================================================================
# 8.7.7  Additional sums of products of two CG               (eq. 205-208)
# ===========================================================================
def r205():
    # Morgan: sum_{l'=0}^{l} (C_{l0,(l'+J)0}^{(l-l'+J)0})^2
    #         / [(2l'-1)(2l-2l'+2J+1)] = -delta_{l0}/(2J+1)
    l = rint(0, 4); J = rint(0, 3)
    tot = S.Zero
    for lp in range(0, int(l) + 1):
        cg = C(l, lp + J, l - lp + J, 0, 0, 0)
        tot += cg ** 2 / ((2 * lp - 1) * (2 * l - 2 * lp + 2 * J + 1))
    rhs = -S.One / (2 * J + 1) if l == 0 else S.Zero
    return tot, rhs


def r206():
    # Morgan: sum_{l'=0}^{l} (1/(2l'+3) - (l+1)/(2l+3) 1/(2l'+1)) /(2l-2l'+2J+1) (CG)^2 = 0
    l = rint(0, 4); J = rint(0, 3)
    tot = S.Zero
    for lp in range(0, int(l) + 1):
        cg = C(l, lp + J, l - lp + J, 0, 0, 0)
        tot += ((S.One / (2 * lp + 3) - (l + 1) / (2 * l + 3) * S.One / (2 * lp + 1))
                / (2 * l - 2 * lp + 2 * J + 1) * cg ** 2)
    return tot, S.Zero


def r207():
    # Din: sum_{i=|c-b|}^{c+b} (2i+1)/(i(i+1)-a(a+1)) (C_{i0,b0}^{c0})^2 = 0
    # (a,b,c non-negative integers, a+b+c odd, |c-b| <= a <= c+b)
    b = rint(0, 4); c = rint(0, 4)
    lo = int(abs(c - b)); hi = int(c + b)
    if lo > hi:
        return None
    a = rint(lo, hi)
    if (a + b + c) % 2 == 0:
        return None
    tot = S.Zero
    for ii in range(lo, hi + 1):
        i = Integer(ii)
        cg = C(i, b, c, 0, 0, 0)
        if cg == 0:                              # sum runs only over non-vanishing CG
            continue                             # (this also skips the pole i=a)
        tot += (2 * i + 1) / (i * (i + 1) - a * (a + 1)) * cg ** 2
    return tot, S.Zero


def r208():
    # Dunlop-Judd: sum_m C_{a m,k 0}^{a m} C_{c M-m,k 0}^{c M-m}
    #  = (-1)^k Pi(a,c)/(2k+1) sqrt((2a-k)!(2c+k+1)!/((2c-k)!(2a+k+1)!)),  a-c >= |M|
    # NOTE: Chap8.tex eq (8.7.208) prints the prefactor as (2a+1)(2c+1); that is
    # an OCR error -- it must be [(2a+1)(2c+1)]^{1/2} = \Fact{a c}, verified here
    # (the printed form is off by exactly sqrt((2a+1)(2c+1))).
    a = rint(1, 4); c = rint(0, int(a))
    k = rint(0, int(2 * c))                       # k integer, <= 2c (and <= 2a)
    if a - c < 0:
        return None
    M = rint(-int(a - c), int(a - c))
    lhs = sum((C(a, k, a, m, 0, m) * C(c, k, c, M - m, 0, M - m)
               for m in proj(a) if abs(M - m) <= c), S.Zero)
    rhs = (ph(k) * Pi(a, c) / (2 * k + 1)
           * sqrt(fac(2 * a - k) * fac(2 * c + k + 1) / (fac(2 * c - k) * fac(2 * a + k + 1))))
    return lhs, rhs


SEC_877 = [
    ("eq 8.7.205  Morgan", r205),
    ("eq 8.7.206  Morgan", r206),
    ("eq 8.7.207  Din", r207),
    ("eq 8.7.208  Dunlop-Judd", r208),
]


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
IMPLEMENTED = [
    ("8.7.1  Sums involving one CG", SEC_871),
    ("8.7.2  Sums of products of two CG", SEC_872),
    ("8.7.3  three CG (one 6j)", SEC_873),
    ("8.7.5  CG and one 6j", SEC_875),
    ("8.7.7  Additional sums of two CG", SEC_877),
]

TODO = [
    ("8.7.4  four CG (one 9j)      eq. 181-191", "some indices still OCR-damaged"),
    ("8.7.6  CG and one 9j         eq. 200-204", "some indices still OCR-damaged"),
]


def run(n, seed):
    random.seed(seed)
    print(f"Section 8.7 sum-rule check -- seed={seed}, up to {n} instances per identity\n")
    all_ok = True
    for title, table in IMPLEMENTED:
        print(f"=== {title} ===")
        for label, fn in table:
            got, bad = 0, None
            for _ in range(n * 400):
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
                print(f"  [SKIP] {label:38s} no valid draws")
                continue
            ok = bad is None
            all_ok &= ok
            print(f"  [{'OK  ' if ok else 'FAIL'}] {label:38s} {got} instances")
            if bad:
                print(f"         lhs={bad[0]}   rhs={bad[1]}")
        print()

    print("Not yet implemented (OCR-damaged; need scan reconstruction):")
    for title, note in TODO:
        print(f"    {title:44s} {note}")
    print()
    print("ALL IMPLEMENTED IDENTITIES HOLD" if all_ok else "SOME IDENTITIES FAILED -- see above")
    return all_ok


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Check Section 8.7 sum rules.")
    p.add_argument("--n", type=int, default=8, help="instances per identity")
    p.add_argument("--seed", type=int, default=20260805, help="RNG seed")
    args = p.parse_args()
    raise SystemExit(0 if run(args.n, args.seed) else 1)
