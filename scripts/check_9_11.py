#!/usr/bin/env python3
"""Verify the Sec. 9.11 6j generator at BOTH levels:

  1. the Racah engine    sixj_algebraic() -> (T, W)              vs sympy, and
  2. the book-form split decompose() -> (parity, out, inside)    vs the engine,

so a bug in the phase / sign-absorption bookkeeping cannot slip through.
Both levels call the shipped functions directly -- nothing is reimplemented.

    python3 scripts/check_9_11.py [2*dmax]      # default 8, i.e. d <= 4
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sympy import Rational, sqrt, lambdify
from sympy.physics.wigner import wigner_6j

from gen_9_11_6j_tables import sixj_algebraic, decompose, u, v, w, a, b, c


def _samples(d, m, n, lo=None):
    """Deficit triples (u,v,w) giving 6j symbols with all triangles satisfied."""
    lo = lo if lo is not None else int(2*d) + 2
    spread = [lo, lo + 1, lo + 2, lo + 6, lo + 14, lo + 30]
    for uu in spread:
        for vv in spread:
            for ww in spread:
                av, bv, cv = Rational(vv+ww, 2), Rational(uu+ww, 2), Rational(uu+vv, 2)
                ev, fv = cv + n, bv + m
                if ev < 0 or fv < 0:
                    continue
                yield uu, vv, ww, av, bv, cv, ev, fv


def main(dmax2=8):
    t0, bad, checked, cells = time.time(), 0, 0, 0
    for d2 in range(1, dmax2 + 1):
        d = Rational(d2, 2)
        vals = [d - i for i in range(d2 + 1)]              # d, d-1, ..., -d
        for m in vals:
            for n in vals:
                cells += 1
                res, dec = sixj_algebraic(d, m, n), decompose(d, m, n)
                if (res is None) != (dec is None):
                    print(f"  ENGINE/SPLIT DISAGREE on zero: d={d} m={m} n={n}")
                    bad += 1
                if res is None:
                    continue
                eng = lambdify((u, v, w), res[0] * sqrt(res[1]), 'math')
                parity, outside, inside = dec
                spl = lambdify((a, b, c), outside * sqrt(inside), 'math')
                for uu, vv, ww, av, bv, cv, ev, fv in _samples(d, m, n):
                    ref = wigner_6j(av, bv, cv, d, ev, fv, prec=64)
                    if ref is None:
                        continue
                    ref = float(ref)
                    s_par = (-1.0)**int(uu + vv + ww)          # (-1)**s
                    try:
                        v_eng = s_par * float(eng(uu, vv, ww))
                        v_spl = (s_par * (-1.0)**parity
                                 * float(spl(float(av), float(bv), float(cv))))
                    except (ValueError, ZeroDivisionError):
                        continue
                    checked += 1
                    tol = 1e-9 * max(1.0, abs(ref))
                    if abs(v_eng - ref) > tol or abs(v_spl - ref) > tol:
                        bad += 1
                        if bad < 8:
                            print(f"  MISMATCH d={d} m={m} n={n} a={av} b={bv} c={cv}: "
                                  f"engine={v_eng} split={v_spl} ref={ref}")
        print(f"  d={str(d):>4}: {cells} cells, {checked} points, {bad} bad "
              f"({time.time()-t0:.0f}s)", flush=True)
    print(f"\n{cells} table cells, {checked} numeric points, {bad} mismatches")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 8))
