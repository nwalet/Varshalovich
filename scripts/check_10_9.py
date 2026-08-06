#!/usr/bin/env python3
r"""
Checks for Section 10.9 (explicit 9j values for special argument values) of
Chapter 10, Varshalovich, Moskalev & Khersonskii.  Built incrementally.

  eq 10.9.1   one argument = 0  -> 6j
  eq 10.9.2   permuted forms all equal (chain)
  eq 10.9.3   two arguments = 0
  eq 10.9.4   {a b c; d e f; 0 0 0}
  eq 10.9.5   {0 b c; d 0 f; g h 0}
  eq 10.9.6   one triad = (g g 1) -> 6j

Usage:  python3 check_10_9.py
"""
import math
from sympy import Rational, S, sqrt
from sympy.physics.wigner import wigner_6j, wigner_9j

H = Rational(1, 2)


def tri(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) == int(a + b + c)


def v9(t):
    a, b, c, d, e, f, g, h, j = t
    return all([tri(a, b, c), tri(d, e, f), tri(g, h, j),
                tri(a, d, g), tri(b, e, h), tri(c, f, j)])


def w9(*t):
    return wigner_9j(*t) if v9(t) else S.Zero


def w6(a, b, c, d, e, f):
    if not (tri(a, b, c) and tri(a, e, f) and tri(d, b, f) and tri(d, e, c)):
        return S.Zero
    return wigner_6j(a, b, c, d, e, f)


def close(u, v):
    d = complex((S(u) - S(v)).evalf(30))
    return math.isfinite(d.real) and abs(d) < 1e-13


R = [0, H, 1, Rational(3, 2), 2]
Rp = [H, 1, Rational(3, 2), 2]


def run():
    print("Section 10.9 checks (zeros + unity)\n")
    ok = True

    # eq 10.9.1 : {a b c; d e f; g h 0}, nonzero c=f, g=h
    g1 = b1 = 0
    for a in Rp:
        for b in Rp:
            for c in Rp:
                for d in Rp:
                    for e in Rp:
                        for g in Rp:
                            t = (a, b, c, d, e, c, g, g, 0)
                            if not v9(t):
                                continue
                            lhs = w9(*t)
                            rhs = (-1) ** (b + c + d + g) / sqrt((2 * c + 1) * (2 * g + 1)) * w6(a, b, c, e, d, g)
                            g1 += 1
                            b1 += 0 if close(lhs, rhs) else 1
    print(f"  [{'OK  ' if g1 and not b1 else 'FAIL'}] eq 10.9.1  {{.. ; g h 0}} -> 6j   ({g1 - b1}/{g1})")
    ok &= g1 > 0 and b1 == 0

    # eq 10.9.4 : {a b c; d e f; 0 0 0}, nonzero a=d,b=e,c=f
    g4 = b4 = 0
    for a in Rp:
        for b in Rp:
            for c in Rp:
                t = (a, b, c, a, b, c, 0, 0, 0)
                if not v9(t):
                    continue
                rhs = 1 / sqrt((2 * a + 1) * (2 * b + 1) * (2 * c + 1))
                g4 += 1
                b4 += 0 if close(w9(*t), rhs) else 1
    print(f"  [{'OK  ' if g4 and not b4 else 'FAIL'}] eq 10.9.4  {{.. ; 0 0 0}}          ({g4 - b4}/{g4})")
    ok &= g4 > 0 and b4 == 0

    # eq 10.9.3 : {a b c; d 0 f; g h 0}; nonzero d=f, b=h, c=f, g=h => b=c=d=f=g=h
    g3 = b3 = 0
    for a in Rp:
        for b in Rp:
            t = (a, b, b, b, 0, b, b, b, 0)
            if not v9(t):
                continue
            rhs = (-1) ** (a - b - b) / ((2 * b + 1) * (2 * b + 1))
            g3 += 1
            b3 += 0 if close(w9(*t), rhs) else 1
    print(f"  [{'OK  ' if g3 and not b3 else 'FAIL'}] eq 10.9.3  two zeros            ({g3 - b3}/{g3})")
    ok &= g3 > 0 and b3 == 0

    # eq 10.9.5 : {0 b c; d 0 f; g h 0}; nonzero b=c=d=f=g=h
    g5 = b5 = 0
    for b in Rp:
        t = (0, b, b, b, 0, b, b, b, 0)
        if not v9(t):
            continue
        rhs = (-1) ** (2 * b) / (2 * b + 1) ** 2
        g5 += 1
        b5 += 0 if close(w9(*t), rhs) else 1
    print(f"  [{'OK  ' if g5 and not b5 else 'FAIL'}] eq 10.9.5  {{0 b c; d 0 f; g h 0}} ({g5 - b5}/{g5})")
    ok &= g5 > 0 and b5 == 0

    # eq 10.9.6 : {a b c; d e c; g g 1} -> 6j
    g6 = b6 = 0
    for a in Rp:
        for b in Rp:
            for c in Rp:
                for d in Rp:
                    for e in Rp:
                        for g in Rp:
                            if g < H:
                                continue
                            t = (a, b, c, d, e, c, g, g, 1)
                            if not v9(t):
                                continue
                            denom = (2 * g + 2) * (2 * g + 1) * (2 * g) * (2 * c + 2) * (2 * c + 1) * (2 * c)
                            if denom == 0:
                                continue
                            rhs = ((-1) ** (b + d + g + c) * 2
                                   * ((a - d) * (a + d + 1) - (b - e) * (b + e + 1)) / sqrt(denom)
                                   * w6(a, b, c, e, d, g))
                            g6 += 1
                            b6 += 0 if close(w9(*t), rhs) else 1
    print(f"  [{'OK  ' if g6 and not b6 else 'FAIL'}] eq 10.9.6  {{.. ; g g 1}} -> 6j    ({g6 - b6}/{g6})")
    ok &= g6 > 0 and b6 == 0

    # eq 10.9.2 : the 8 permuted 9j are all equal
    g2 = b2 = 0
    for c in Rp:
        for e in Rp:
            for b in Rp:
                for g in Rp:
                    for d in Rp:
                        for a in Rp:
                            ref = w9(0, c, c, g, e, b, g, d, a)
                            if not v9((0, c, c, g, e, b, g, d, a)):
                                continue
                            forms = [(c, 0, c, d, g, a, e, g, b), (g, g, 0, e, d, c, b, a, c),
                                     (g, b, e, 0, c, c, g, a, d), (a, g, d, c, 0, c, b, g, e),
                                     (b, a, c, g, g, 0, e, d, c), (c, e, d, c, b, a, 0, g, g),
                                     (d, c, e, a, c, b, g, 0, g)]
                            for fm in forms:
                                g2 += 1
                                b2 += 0 if close(w9(*fm), ref) else 1
    print(f"  [{'OK  ' if g2 and not b2 else 'FAIL'}] eq 10.9.2  8 permuted 9j equal  ({g2 - b2}/{g2})")
    ok &= g2 > 0 and b2 == 0

    # eq 10.9.7 : {a b c; d e c; g+1 g 1} recursion into two 6j
    g7 = b7 = 0
    for a in Rp:
        for b in Rp:
            for c in Rp:
                for d in Rp:
                    for e in Rp:
                        for g in Rp:
                            t = (a, b, c, d, e, c, g + 1, g, 1)
                            if not v9(t):
                                continue
                            fac = (-1) ** (b + d + g + c) * sqrt(Rational((2 * g + 3) * (2 * g + 2) * (2 * g + 1) * (2 * c + 2) * (2 * c + 1) * (2 * c), 2))
                            lhs = w9(*t) * fac

                            def rad(*xs):
                                p = S.One
                                for x in xs:
                                    if x < 0:
                                        return None
                                    p *= x
                                return sqrt(p)
                            r1 = rad(b - e + g + 1, -b + e + g + 1, b + e + g + 2, b + e - g)
                            r2 = rad(a - d + g + 1, -a + d + g + 1, a + d + g + 2, a + d - g)
                            if r1 is None or r2 is None:
                                continue
                            rhs = r1 * w6(a, d, g + 1, e, b, c) + r2 * w6(a, d, g, e, b, c)
                            g7 += 1
                            b7 += 0 if close(lhs, rhs) else 1
    print(f"  [{'OK  ' if g7 and not b7 else 'FAIL'}] eq 10.9.7  {{.. ; g+1 g 1}} recursion ({g7 - b7}/{g7})")
    ok &= g7 > 0 and b7 == 0

    print("\nALL CHECKED 10.9 FORMS PASS" if ok else "\nSOME 10.9 CHECKS FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
