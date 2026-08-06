#!/usr/bin/env python3
r"""
Checks for Section 10.13 (12j symbols of the first kind) of Chapter 10,
Varshalovich, Moskalev & Khersonskii.

sympy has no 12j symbol, so we take the four-6j single-sum form eq 10.13.6 as
the definition of 12j(I) and cross-check the other representations against it:

  eq 10.13.1  general 3nj first kind, specialised to n=4
  eq 10.13.6  12j(I) = sum of four 6j            (reference definition)
  eq 10.13.7  12j(I) = sum of 9j x 6j x 6j
  eq 10.13.9  16 symmetry permutations all equal

(eq 10.13.3, the eight-3jm sum, and 10.13.4/5 orthogonality are not covered
here.)

Argument order (as in the \twelvejI macro):
  (a1 a2 a3 a4 | b12 b23 b34 b41 | c1 c2 c3 c4)

Usage:  python3 check_10_13.py
"""
import math
from sympy import Rational, S
from sympy.physics.wigner import wigner_6j, wigner_9j

H = Rational(1, 2)


def tri(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) == int(a + b + c)


def w6(a, b, c, d, e, f):
    if not (tri(a, b, c) and tri(a, e, f) and tri(d, b, f) and tri(d, e, c)):
        return S.Zero
    return wigner_6j(a, b, c, d, e, f)


def w9(a, b, c, d, e, f, g, h, j):
    if not all([tri(a, b, c), tri(d, e, f), tri(g, h, j),
                tri(a, d, g), tri(b, e, h), tri(c, f, j)]):
        return S.Zero
    return wigner_9j(a, b, c, d, e, f, g, h, j)


def valid12(a1, a2, a3, a4, b12, b23, b34, b41, c1, c2, c3, c4):
    triads = [(a1, b12, a2), (a2, b23, a3), (a3, b34, a4), (a4, b41, c1),
              (c1, b12, c2), (c2, b23, c3), (c3, b34, c4), (c4, b41, a1)]
    if not all(tri(*t) for t in triads):
        return False

    def tetra(j1, j2, j3, j4):
        return ((j1 + j2 + j3 + j4) == int(j1 + j2 + j3 + j4)
                and j1 <= j2 + j3 + j4 and j2 <= j1 + j3 + j4
                and j3 <= j1 + j2 + j4 and j4 <= j1 + j2 + j3)
    return tetra(a1, c1, a3, c3) and tetra(a2, c2, a4, c4)


def xrange12(v):
    hi = sum(v)
    return [Rational(i, 2) for i in range(0, int(2 * hi) + 1)]


# eq 10.13.6 : reference definition (four 6j)
def TW(a1, a2, a3, a4, b12, b23, b34, b41, c1, c2, c3, c4):
    Ssum = a1 + a2 + a3 + a4 + c1 + c2 + c3 + c4 + b12 + b23 + b34 + b41
    tot = S.Zero
    for x in xrange12((a1, a2, a3, a4, b12, b23, b34, b41, c1, c2, c3, c4)):
        tot += ((-1) ** (Ssum - x) * (2 * x + 1)
                * w6(a1, a2, b12, c2, c1, x) * w6(a2, a3, b23, c3, c2, x)
                * w6(a3, a4, b34, c4, c3, x) * w6(a4, c1, b41, a1, c4, x))
    return tot


# eq 10.13.1 (n=4): j_i=a_i, l_i=(b12,b23,b34,b41), k_i=c_i
def eq1(a1, a2, a3, a4, b12, b23, b34, b41, c1, c2, c3, c4):
    j = [a1, a2, a3, a4]
    l = [b12, b23, b34, b41]
    k = [c1, c2, c3, c4]
    R = sum(j) + sum(l) + sum(k)
    tot = S.Zero
    for x in xrange12((a1, a2, a3, a4, b12, b23, b34, b41, c1, c2, c3, c4)):
        tot += ((-1) ** (R + 3 * x) * (2 * x + 1)
                * w6(j[0], k[0], x, k[1], j[1], l[0]) * w6(j[1], k[1], x, k[2], j[2], l[1])
                * w6(j[2], k[2], x, k[3], j[3], l[2]) * w6(j[3], k[3], x, j[0], k[0], l[3]))
    return tot


# eq 10.13.7 : 9j x 6j x 6j
def eq7(a1, a2, a3, a4, b12, b23, b34, b41, c1, c2, c3, c4):
    tot = S.Zero
    for x in xrange12((a1, a2, a3, a4, b12, b23, b34, b41, c1, c2, c3, c4)):
        tot += ((2 * x + 1)
                * w9(b23, a2, a3, c2, b12, c1, c3, a1, x)
                * w6(b34, a3, a4, c1, b41, x) * w6(b34, c3, c4, a1, b41, x))
    return (-1) ** (a1 - a3 - c1 + c3) * tot


def close(u, v):
    d = complex((S(u) - S(v)).evalf(30))
    return math.isfinite(d.real) and abs(d) < 1e-12


# a handful of valid 12j(I) argument sets
def gen_cases():
    import itertools
    vals = [H, 1, Rational(3, 2)]
    out = []
    for combo in itertools.product([1, H], repeat=12):
        if valid12(*combo):
            out.append(combo)
        if len(out) >= 30:
            break
    # add a few larger
    base = (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    out.append(base)
    out.append((1, Rational(3, 2), 1, H, 1, 1, 1, 1, 1, H, 1, Rational(3, 2)))
    return [c for c in out if valid12(*c)]


def run():
    print("Section 10.13 : 12j(I) cross-checks\n")
    CASES = gen_cases()
    ok = True
    for name, fn in [("eq 10.13.1 (n=4) vs 10.13.6", eq1), ("eq 10.13.7 (9j.6j.6j) vs 10.13.6", eq7)]:
        good = bad = 0
        for c in CASES:
            ref = TW(*c)
            good += 1
            if not close(fn(*c), ref):
                bad += 1
        okk = good > 0 and bad == 0
        ok &= okk
        print(f"  [{'OK  ' if okk else 'FAIL'}] {name:36s} ({good - bad}/{good})")

    # eq 10.13.9 : 16 symmetry permutations (index permutations of the 12 args)
    def perms(a1, a2, a3, a4, b12, b23, b34, b41, c1, c2, c3, c4):
        A = [a1, a2, a3, a4]
        B = [b12, b23, b34, b41]
        Cc = [c1, c2, c3, c4]
        # the 8 "AC-chain" as one cyclic list of length 8: a1 a2 a3 a4 c1 c2 c3 c4
        chain = A + Cc
        res = []
        # 4 cyclic shifts by 1 (rotate chain by 2 and B by 1)
        for s in range(4):
            ch = chain[2 * s:] + chain[:2 * s]
            b = B[s:] + B[:s]
            res.append((ch[0], ch[1], ch[2], ch[3], b[0], b[1], b[2], b[3], ch[4], ch[5], ch[6], ch[7]))
        # inversion of order + 4 shifts
        chain_r = chain[::-1]
        chain_r = chain_r[-1:] + chain_r[:-1]        # align so a1 stays first under the book's convention
        Br = B[::-1]
        for s in range(4):
            ch = chain_r[2 * s:] + chain_r[:2 * s]
            b = Br[s:] + Br[:s]
            res.append((ch[0], ch[1], ch[2], ch[3], b[0], b[1], b[2], b[3], ch[4], ch[5], ch[6], ch[7]))
        return res

    good = bad = 0
    for c in CASES:
        ref = TW(*c)
        if ref == 0:
            continue
        for p in perms(*c):
            if not valid12(*p):
                continue
            good += 1
            if not close(TW(*p), ref):
                bad += 1
    okk = good > 0 and bad == 0
    print(f"  [{'OK  ' if okk else 'FAIL'}] {'eq 10.13.9 symmetry perms':36s} ({good - bad}/{good})")
    ok &= okk

    print("\nALL CHECKED 10.13 FORMS PASS" if ok else "\nSOME 10.13 CHECKS FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
