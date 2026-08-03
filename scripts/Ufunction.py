#!/usr/bin/env python3
"""U-function U^J_{M M'}(omega; Theta, Phi) in algebraic form -- VMK convention.

U^J_{M M'}(omega; Theta, Phi) = <J M| exp(-i omega n.J) |J M'> is the matrix of a
rotation by the angle ``omega`` about the axis

    n(Theta, Phi) = (sin Theta cos Phi, sin Theta sin Phi, cos Theta),

in the convention of Varshalovich-Moskalev-Khersonskii (VMK), Sec. 4.5.  It
reduces to the small Wigner d-function for a rotation about the y-axis:

    U^J_{M M'}(omega; pi/2, pi/2) == d^J_{M M'}(omega)       (see wigner_d.py).

The closed form has a REMOVABLE singularity at Theta = 0 (there v = 0 and the
1/v**2 terms cancel analytically); use the symbolic result, or keep Theta != 0
for purely numerical evaluation.  Pass half-integer J, M, M' as ``Rational(1, 2)``
etc. so the result stays algebraic.

This is the SymPy analogue of ``Ufunction.wl``.  Verified against the matrix
exponential exp(-i omega n.J) (via numpy/scipy) for j = 1/2 .. 2, and against
Rotation.d for the y-rotation slice.  Run ``--selftest`` to redo the checks.

Requires: sympy (numpy + scipy optional, only for the full self-test).  Usage:
    python3 Ufunction.py               # demo
    python3 Ufunction.py --selftest    # re-validate
"""
from sympy import (symbols, I, sqrt, factorial, exp, cos, sin, pi, simplify,
                   nsimplify, Matrix, Rational)
from sympy.physics.quantum.spin import Rotation

omega, Theta, Phi = symbols('omega Theta Phi', real=True)


def u_function(j, m, mp, w=omega, th=Theta, ph=Phi):
    """Algebraic U^J_{m,mp}(w; th, ph) in the VMK convention."""
    j, m, mp = nsimplify(j), nsimplify(m), nsimplify(mp)
    ampp = abs(m + mp)
    if m + mp >= 0:
        uu = cos(w/2) - I*sin(w/2)*cos(th)
        maxs = min(j - m, j - mp)
        term = lambda s: 1/(factorial(s)*factorial(s + ampp)
                            * factorial(j - m - s)*factorial(j - mp - s))
    else:
        uu = cos(w/2) + I*sin(w/2)*cos(th)
        maxs = min(j + m, j + mp)
        term = lambda s: 1/(factorial(s)*factorial(s + ampp)
                            * factorial(j + m - s)*factorial(j + mp - s))
    v = sin(w/2)*sin(th)
    pref = ((-I*v)**(2*j - ampp) * uu**ampp * exp(-I*(m - mp)*ph)
            * sqrt(factorial(j + m)*factorial(j - m)
                   * factorial(j + mp)*factorial(j - mp)))
    tot = sum(term(s)*(1 - v**-2)**s for s in range(int(maxs) + 1))
    return simplify(pref*tot)


def u_function_matrix(j, w=omega, th=Theta, ph=Phi):
    """Full U^J(w; th, ph) as a SymPy Matrix; rows M and cols M' run J..-J."""
    j = nsimplify(j)
    order = [j - k for k in range(int(2*j) + 1)]
    return Matrix([[u_function(j, m, mp, w, th, ph) for mp in order]
                   for m in order])


def _selftest():
    js = [Rational(1, 2), 1, Rational(3, 2), 2]
    worst = 0.0
    try:
        import numpy as np
        from scipy.linalg import expm

        def spin_mats(j):
            ms = [j - k for k in range(int(2*j) + 1)]
            n = len(ms)
            Jp = np.zeros((n, n), complex)
            for k in range(n):
                if k - 1 >= 0:
                    Jp[k - 1, k] = np.sqrt(j*(j + 1) - ms[k]*(ms[k] + 1))
            return ((Jp + Jp.conj().T)/2, (Jp - Jp.conj().T)/(2j),
                    np.diag(ms).astype(complex), ms)

        def ref(j, m, mp, w, th, ph):
            Jx, Jy, Jz, ms = spin_mats(j)
            nn = [np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)]
            U = expm(-1j*w*(nn[0]*Jx + nn[1]*Jy + nn[2]*Jz))
            return U[ms.index(m), ms.index(mp)]

        pts = [(0.7, 1.1, 0.3), (2.0, 0.4, 1.7), (1.3, 2.5, -0.9)]
        for j in js:
            order = [j - k for k in range(int(2*j) + 1)]
            for m in order:
                for mp in order:
                    e = u_function(j, m, mp)
                    for (w, th, ph) in pts:
                        val = complex(e.subs({omega: w, Theta: th, Phi: ph}).evalf(30))
                        worst = max(worst, abs(val - ref(float(j), float(m),
                                                         float(mp), w, th, ph)))
        print(f"self-test vs scipy expm(-i w n.J), j=1/2..2: "
              f"worst |diff| = {worst:.2e}")
    except ImportError:
        print("self-test: numpy/scipy not available, skipping matrix-exponential check")

    # analytic cross-check (sympy only): U(w; pi/2, pi/2) == d^J(w)
    beta = symbols('beta', real=True)
    ok = True
    for j in js:
        order = [j - k for k in range(int(2*j) + 1)]
        for m in order:
            for mp in order:
                du = u_function(j, m, mp).subs({Theta: pi/2, Phi: pi/2})
                dd = Rotation.d(j, m, mp, omega).doit()
                for w in (0.4, 1.3, 2.7):
                    if abs(complex(du.subs(omega, w))
                           - complex(dd.subs(beta, w).subs(omega, w))) > 1e-12:
                        ok = False
    print("cross-check U(w; pi/2, pi/2) == d^J(w):", "OK" if ok else "FAIL")
    return (worst < 1e-10) and ok


if __name__ == '__main__':
    import sys
    if '--selftest' in sys.argv:
        raise SystemExit(0 if _selftest() else 1)

    # --- specific quantum numbers ---
    print("U^{1/2}_{1/2,1/2}(w; Th, Ph) =",
          u_function(Rational(1, 2), Rational(1, 2), Rational(1, 2)))
    print("U^{1}_{1,0}(w; Th, Ph)       =", u_function(1, 1, 0))
    print()
    print("U(w; pi/2, pi/2) for j=1  (== d^1(w)):")
    print(u_function_matrix(1).subs({Theta: pi/2, Phi: pi/2}).applyfunc(simplify))
