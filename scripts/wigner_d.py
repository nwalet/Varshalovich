#!/usr/bin/env python3
"""Wigner (small) d-function d^J_{M M'}(beta) in algebraic form -- VMK convention.

Uses SymPy's ``sympy.physics.quantum.spin.Rotation.d``, whose convention matches
Varshalovich-Moskalev-Khersonskii (VMK):

    d^{1/2}_{1/2,-1/2}(beta) = -sin(beta/2)

and the d^1 matrix reproduces the VMK table.  This has been checked against the
explicit Wigner sum formula for j = 1/2, 1, 3/2, 2, 5/2 and all m, m'
(proof-by-points at 15 angles, 30-digit precision).  Run ``--selftest`` to redo it.

Requires: sympy.  Usage:
    python3 wigner_d.py               # demo
    python3 wigner_d.py --selftest    # re-validate against the Wigner formula
"""
from sympy import (symbols, Rational, simplify, Matrix, nsimplify, cos, sin,
                   sqrt, factorial, S, pi)
from sympy.physics.quantum.spin import Rotation

beta = symbols('beta', real=True)


def wigner_d(j, m, mp, angle=beta):
    """Algebraic d^j_{m,mp}(angle) in VMK / Wigner convention.

    j, m, mp may be integers or half-integers; pass half-integers as
    ``Rational(1, 2)`` (or the string ``'1/2'``) so the result stays symbolic.
    """
    j, m, mp = nsimplify(j), nsimplify(m), nsimplify(mp)
    return simplify(Rotation.d(j, m, mp, angle).doit())


def wigner_d_matrix(j, angle=beta):
    """Full d^j(angle) as a SymPy Matrix; rows m and cols m' run j, j-1, ..., -j."""
    j = nsimplify(j)
    order = [j - k for k in range(int(2 * j) + 1)]
    return Matrix([[wigner_d(j, m, mp, angle) for mp in order] for m in order])


def _reference_d(j, m, mp):
    """Independent Wigner explicit sum formula (VMK eq. 4.3.1), for validation."""
    j, m, mp = nsimplify(j), nsimplify(m), nsimplify(mp)
    pref = sqrt(factorial(j + m) * factorial(j - m) *
                factorial(j + mp) * factorial(j - mp))
    tot = S.Zero
    for k in range(int(max(0, m - mp)), int(min(j + m, j - mp)) + 1):
        tot += ((-1)**k * cos(beta / 2)**(2*j - 2*k + m - mp)
                * sin(beta / 2)**(2*k - m + mp)
                / (factorial(k) * factorial(j + m - k)
                   * factorial(j - mp - k) * factorial(k - m + mp)))
    return pref * tot


def _selftest():
    """d^j is a degree-2j trig polynomial in beta/2, so vanishing of the
    difference at > 2*(2j)+1 distinct points proves the identity."""
    pts = [Rational(p, 17) * pi for p in range(1, 16)]      # 15 distinct angles
    worst = 0.0
    for j2 in range(1, 6):                                   # 2j = 1..5
        j = Rational(j2, 2)
        order = [j - k for k in range(int(2 * j) + 1)]
        for m in order:
            for mp in order:
                e = Rotation.d(j, m, mp, beta).doit() - _reference_d(j, m, mp)
                for p in pts:
                    worst = max(worst, abs(complex(e.subs(beta, p).evalf(30))))
    ok = worst < 1e-25
    print(f"self-test j=1/2..5/2 vs Wigner formula: "
          f"{'IDENTICAL' if ok else 'FAILED'}  (worst |diff| = {worst:.2e})")
    return ok


if __name__ == '__main__':
    import sys
    if '--selftest' in sys.argv:
        raise SystemExit(0 if _selftest() else 1)

    # --- specific angular momentum / projection quantum numbers ---
    print("d^{1/2}_{1/2,-1/2}(beta) =",
          wigner_d(Rational(1, 2), Rational(1, 2), Rational(-1, 2)))
    print("d^{2}_{1,0}(beta)        =", wigner_d(2, 1, 0))
    print()
    print("d^1(beta) =")
    print(wigner_d_matrix(1))
