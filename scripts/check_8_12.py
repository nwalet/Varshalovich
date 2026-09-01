#!/usr/bin/env python3
"""Verify the Sec. 8.12 CG generator at BOTH levels:

  1. the Racah engine   cg_algebraic() -> (T, W)          against sympy's CG, and
  2. the book-form split decompose() -> (sign, out, in)   against the engine,

so a bug in the sign/absorption bookkeeping cannot slip through.  Both levels
call the shipped functions directly -- nothing here reimplements them.

    python3 scripts/check_8_12.py [2*bmax]      # default 10, i.e. b <= 5
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sympy import Rational, sqrt, lambdify
from sympy.physics.quantum.cg import CG

from gen_8_12_cg_tables import cg_algebraic, decompose, p, q, c, g


def main(bmax2=10):
    t0, bad, checked, cells, nonzero = time.time(), 0, 0, 0, 0
    for b2 in range(1, bmax2 + 1):
        b = Rational(b2, 2)
        vals = [b - i for i in range(b2 + 1)]              # b, b-1, ..., -b
        for k in vals:
            for beta in vals:
                cells += 1
                res, dec = cg_algebraic(b, beta, k), decompose(b, beta, k)
                if (res is None) != (dec is None):
                    print(f"  ENGINE/SPLIT DISAGREE on zero: b={b} beta={beta} k={k}")
                    bad += 1
                if res is None:
                    continue
                nonzero += 1
                eng = lambdify((p, q), res[0] * sqrt(res[1]), 'math')
                sign, outside, inside = dec
                ren = lambdify((c, g), sign * outside * sqrt(inside), 'math')
                for a2 in range(b2 + 2, b2 + 8):
                    a = Rational(a2, 2)
                    for al2 in range(-a2, a2 + 1, 2):
                        al = Rational(al2, 2)
                        cc, gg = a + k, al + beta
                        if cc <= 0 or abs(gg) > cc or abs(al) > a:
                            continue
                        ref = float(CG(a, al, b, beta, cc, gg).doit())
                        try:
                            v_eng = float(eng(float(a + al), float(a - al)))
                            v_ren = float(ren(float(cc), float(gg)))
                        except (ValueError, ZeroDivisionError):
                            continue
                        checked += 1
                        tol = 1e-9 * max(1.0, abs(ref))
                        if abs(v_eng - ref) > tol or abs(v_ren - ref) > tol:
                            bad += 1
                            if bad < 8:
                                print(f"  MISMATCH b={b} beta={beta} k={k} a={a} "
                                      f"al={al}: engine={v_eng} split={v_ren} ref={ref}")
        print(f"  b={str(b):>4}: {cells} cells ({nonzero} non-zero), "
              f"{checked} points, {bad} bad ({time.time()-t0:.0f}s)", flush=True)
    print(f"\n{cells} table cells, {nonzero} non-zero, "
          f"{checked} numeric points, {bad} mismatches")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 10))
