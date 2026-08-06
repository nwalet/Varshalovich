#!/usr/bin/env python3
r"""
Checks for Section 9.8 (sums involving 6j symbols) of Chapter 9, Varshalovich,
Moskalev & Khersonskii.

  eq 9.8.1  sum_X (2X+1){a b X; a b c} = (-1)^{2c}{abc}
  eq 9.8.2  sum_X (-1)^{a+b+X}(2X+1){a b X; b a c} = delta_{c0} sqrt((2a+1)(2b+1))
  eq 9.8.3  sum_X (2X+1){a b X;c d p}{a b X;c d q} = delta_pq {adp}{bcp}/(2p+1)
  eq 9.8.4  sum_X (-1)^{p+q+X}(2X+1){a b X;c d p}{a b X;d c q} = {a c q; b d p}
  eq 9.8.5  sum_X (-1)^{2X}(2X+1){a b X;c d p}{c d X;e f q}{e f X;a b r} = 9j
  eq 9.8.6  sum_X (-1)^{R+X}(2X+1){a b X;c d p}{c d X;e f q}{e f X;b a r}
              = {p q r;e a d}{p q r;f b c}

Not covered: 9.8.7/9.8.8 (12j symbols, not in sympy).

Usage:  python3 check_9_8.py
"""
import math
from sympy import Rational, S
from sympy.physics.wigner import wigner_6j, wigner_9j

H = Rational(1, 2)


def tri(a, b, c):
    return (a >= 0 and b >= 0 and c >= 0 and abs(a - b) <= c <= a + b
            and (a + b + c) == int(a + b + c))


def T(a, b, c):
    return 1 if tri(a, b, c) else 0


def v6(a, b, c, d, e, f):
    return tri(a, b, c) and tri(c, d, e) and tri(a, e, f) and tri(b, d, f)


def w6(a, b, c, d, e, f):
    return wigner_6j(a, b, c, d, e, f) if v6(a, b, c, d, e, f) else S.Zero


def w9(a, b, c, d, e, f, g, h, i):
    if not all([tri(a, b, c), tri(d, e, f), tri(g, h, i),
                tri(a, d, g), tri(b, e, h), tri(c, f, i)]):
        return S.Zero
    return wigner_9j(a, b, c, d, e, f, g, h, i)


def Xrange(hi):
    return [Rational(i, 2) for i in range(0, int(2 * hi) + 1)]


def close(u, v):
    d = complex((S(u) - S(v)).evalf(30))
    return math.isfinite(d.real) and abs(d) < 1e-13


VJ = [0, H, 1, Rational(3, 2), 2, Rational(5, 2), 3]
HI = 6


def eq1():
    good = bad = 0
    for a in VJ:
        for b in VJ:
            for c in VJ:
                lhs = sum((2 * X + 1) * w6(a, b, X, a, b, c) for X in Xrange(HI))
                rhs = (-1) ** (2 * c) * T(a, b, c)
                good += 1
                if not close(lhs, rhs):
                    bad += 1
    return good, bad


def eq2():
    good = bad = 0
    for a in VJ:
        for b in VJ:
            for c in VJ:
                lhs = sum((-1) ** (a + b + X) * (2 * X + 1) * w6(a, b, X, b, a, c) for X in Xrange(HI))
                rhs = (1 if c == 0 else 0) * (S(2 * a + 1) * (2 * b + 1)) ** H
                good += 1
                if not close(lhs, rhs):
                    bad += 1
    return good, bad


def eq3():
    good = bad = 0
    for a in [H, 1, Rational(3, 2), 2]:
        for b in [H, 1, Rational(3, 2), 2]:
            for c in [H, 1, Rational(3, 2), 2]:
                for d in [H, 1, Rational(3, 2), 2]:
                    for p in Xrange(HI):
                        for q in Xrange(HI):
                            lhs = sum((2 * X + 1) * w6(a, b, X, c, d, p) * w6(a, b, X, c, d, q)
                                      for X in Xrange(HI))
                            rhs = (1 if p == q else 0) * T(a, d, p) * T(b, c, p) / (2 * p + 1)
                            good += 1
                            if not close(lhs, rhs):
                                bad += 1
    return good, bad


def eq4():
    good = bad = 0
    for a in [H, 1, Rational(3, 2), 2]:
        for b in [H, 1, Rational(3, 2), 2]:
            for c in [H, 1, Rational(3, 2), 2]:
                for d in [H, 1, Rational(3, 2), 2]:
                    for p in Xrange(HI):
                        for q in Xrange(HI):
                            lhs = sum((-1) ** (p + q + X) * (2 * X + 1)
                                      * w6(a, b, X, c, d, p) * w6(a, b, X, d, c, q)
                                      for X in Xrange(HI))
                            rhs = w6(a, c, q, b, d, p)
                            good += 1
                            if not close(lhs, rhs):
                                bad += 1
    return good, bad


def eq5():
    good = bad = 0
    P = [H, 1]
    for a in P:
        for b in P:
            for c in P:
                for d in P:
                    for e in P:
                        for f in P:
                            for p in P:
                                for q in P:
                                    for r in P:
                                        lhs = sum((-1) ** (2 * X) * (2 * X + 1)
                                                  * w6(a, b, X, c, d, p) * w6(c, d, X, e, f, q) * w6(e, f, X, a, b, r)
                                                  for X in Xrange(HI))
                                        rhs = w9(a, f, r, d, q, e, p, c, b)
                                        if lhs == 0 and rhs == 0:
                                            continue
                                        good += 1
                                        if not close(lhs, rhs):
                                            bad += 1
    return good, bad


def eq6():
    good = bad = 0
    P = [H, 1]
    for a in P:
        for b in P:
            for c in P:
                for d in P:
                    for e in P:
                        for f in P:
                            for p in P:
                                for q in P:
                                    for r in P:
                                        R = a + b + c + d + e + f + p + q + r
                                        lhs = sum((-1) ** (R + X) * (2 * X + 1)
                                                  * w6(a, b, X, c, d, p) * w6(c, d, X, e, f, q) * w6(e, f, X, b, a, r)
                                                  for X in Xrange(HI))
                                        rhs = w6(p, q, r, e, a, d) * w6(p, q, r, f, b, c)
                                        if lhs == 0 and rhs == 0:
                                            continue
                                        good += 1
                                        if not close(lhs, rhs):
                                            bad += 1
    return good, bad


def run():
    print("Section 9.8 sum-rule checks\n")
    ok = True
    for name, fn in [("eq 9.8.1", eq1), ("eq 9.8.2", eq2), ("eq 9.8.3", eq3),
                     ("eq 9.8.4", eq4), ("eq 9.8.5 (=9j)", eq5), ("eq 9.8.6", eq6)]:
        good, bad = fn()
        okk = good > 0 and bad == 0
        ok &= okk
        print(f"  [{'OK  ' if okk else 'FAIL'}] {name:16s} ({good} cases, {bad} bad)")
    print("\n  (not checked: 9.8.7, 9.8.8 -- 12j symbols)")
    print("\nALL 9.8 CHECKS PASS" if ok else "\nSOME 9.8 CHECKS FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
