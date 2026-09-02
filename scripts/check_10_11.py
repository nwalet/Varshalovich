#!/usr/bin/env python3
"""Verify the Sec. 10.11 9j generator against sympy's wigner_9j.

Checks the shipped decompose() -- sign, the factor left outside the radical and
the absorbed radicand -- so a misclassified sign-definite factor cannot slip
through.  Sample momenta are deliberately ASYMMETRIC as well as balanced: a
wrongly absorbed factor f is rendered |f| instead of f, which differs only
where f < 0, and factors like Z only go negative when one momentum dominates.

    python3 scripts/check_10_11.py [max_2alpha]     # default 6 (alpha <= 3)
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# lambdify compiles a deeply nested tree; the largest cells overflow the default
sys.setrecursionlimit(20000)

from sympy import Rational, sqrt, lambdify
from sympy.physics.wigner import wigner_9j

from gen_10_11_9j_tables import decompose, a, b, c

# the (alpha, beta, gamma) triples of Tables 10.1-10.12
TABLES = [(Rational(1,2), Rational(1,2), 0), (Rational(1,2), Rational(1,2), 1),
          (1, 1, 0), (1, 1, 1),
          (Rational(3,2), Rational(3,2), 0), (Rational(3,2), Rational(3,2), 1),
          (2, 1, 1), (2, 2, 0), (2, 2, 1),
          (Rational(5,2), Rational(3,2), 1), (3, 2, 1)]

# balanced and strongly asymmetric deficit triples (u, v, w)
DEFICITS = [(6, 6, 6), (8, 7, 9), (12, 4, 5), (4, 12, 5), (5, 4, 12),
            (20, 5, 6), (5, 20, 6), (6, 5, 20), (16, 15, 4), (4, 15, 16)]


def main(amax2=6):
    t0, bad, checked, cells = time.time(), 0, 0, 0
    for (al, be, ga) in TABLES:
        if 2*Rational(al) > amax2:
            continue
        for lam in [al - i for i in range(int(2*al) + 1)]:
            for mu in [be - i for i in range(int(2*be) + 1)]:
                for nu in [ga - i for i in range(int(2*ga) + 1)]:
                    cells += 1
                    try:
                        dec = decompose(al, be, ga, lam, mu, nu)
                    except ValueError as e:      # radicands failed to share a radical
                        print(f"  COMBINE FAIL a={al} b={be} g={ga} "
                              f"lam={lam} mu={mu} nu={nu}: {e}")
                        bad += 1
                        continue
                    if dec is None:
                        f = None
                    else:
                        sign, outside, inside = dec
                        f = lambdify((a, b, c), sign * outside * sqrt(inside), 'math')
                    for (uu, vv, ww) in DEFICITS:
                        av, bv, cv = (Rational(vv+ww, 2), Rational(uu+ww, 2),
                                      Rational(uu+vv, 2))
                        if (av+bv+cv) % 1 or (av+lam+bv+mu+cv+nu) % 1:
                            continue
                        ref = wigner_9j(av+lam, bv+mu, cv+nu, av, bv, cv,
                                        al, be, ga, prec=64)
                        ref = 0.0 if ref is None else float(ref)
                        try:
                            got = 0.0 if f is None else float(
                                f(float(av), float(bv), float(cv)))
                        except (ValueError, ZeroDivisionError):
                            continue
                        checked += 1
                        if abs(got - ref) > 1e-9 * max(1.0, abs(ref)):
                            bad += 1
                            if bad < 8:
                                print(f"  MISMATCH al={al} be={be} ga={ga} lam={lam} "
                                      f"mu={mu} nu={nu} abc={av},{bv},{cv}: "
                                      f"{got} vs {ref}")
        print(f"  ({al},{be},{ga}): {cells} cells, {checked} points, {bad} bad "
              f"({time.time()-t0:.0f}s)", flush=True)
    print(f"\n{cells} table cells, {checked} numeric points, {bad} mismatches")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 6))
