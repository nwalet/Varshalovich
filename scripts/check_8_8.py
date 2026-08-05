#!/usr/bin/env python3
r"""
Symbolic check of the generating-function identities in Section 8.8 of
Varshalovich, Moskalev & Khersonskii, "Quantum Theory of Angular Momentum".

Unlike the random-sampling checkers for Secs. 8.4-8.7, these identities reduce
(for fixed a,b,c and, where relevant, gamma) to polynomial / power-series /
single-variable-function identities, so they are verified symbolically with
sympy: expand both sides and confirm the difference is identically zero (or,
for the hypergeometric / d-function cases, evaluate at several points).

Covered
    8.8  eq 8.8.2   products of binomials, one variable per momentum
    8.8  eq 8.8.3  products of binomials, two variables per momentum
    8.8  eq 8.8.4  exponential generating function
    8.8  eq 8.8.5   hypergeometric generating function (Akim-Levin)
    8.8  eq 8.8.6   Wigner d-function generating function
    8.8  eq 8.8.7   Schwinger generating function for C_{a0,b0}^{c0}

Not covered
    8.8  eq 8.8.1   Regge determinant to power J -- a sum over all R-symbols at
                  fixed J; checkable only by enumerating them, cumbersome.

Usage:  python3 check_8_8.py
"""
from sympy import (symbols, sqrt, Rational, factorial as fac, expand, cos, sin,
                   hyper, hyperexpand, S)
from sympy.physics.wigner import clebsch_gordan as CG
from sympy.physics.quantum.spin import Rotation

t, t1, t2, t3, th = symbols('t t1 t2 t3 th')
u1, u2, v1, v2, w1, w2 = symbols('u1 u2 v1 v2 w1 w2')


def rng(j):
    return [-j + i for i in range(int(2 * j) + 1)]


def _sum_abg(a, b, c, term):
    """sum term(al, be, ga) over al,be with ga = -(al+be), |ga| <= c."""
    tot = S.Zero
    for al in rng(a):
        for be in rng(b):
            ga = -(al + be)
            if abs(ga) > c:
                continue
            tot += term(al, be, ga)
    return tot


def eq210(a, b, c):
    J = a + b + c
    lhs = (t1 - t2) ** (J - 2 * c) * (t2 - t3) ** (J - 2 * a) * (t3 - t1) ** (J - 2 * b)
    pref = sqrt(Rational(fac(J + 1) * fac(J - 2 * a) * fac(J - 2 * b) * fac(J - 2 * c), 2 * c + 1))
    rhs = pref * _sum_abg(a, b, c, lambda al, be, ga:
        S.NegativeOne ** (a - b - ga) * t1 ** (a + al) * t2 ** (b + be) * t3 ** (c + ga)
        / sqrt(fac(a + al) * fac(a - al) * fac(b + be) * fac(b - be) * fac(c + ga) * fac(c - ga))
        * CG(a, b, c, al, be, -ga))
    return expand(lhs - rhs) == 0


def eq210b(a, b, c):
    J = a + b + c
    lhs = ((v1 * u2 - u1 * v2) ** (J - 2 * c) * (w1 * v2 - w2 * v1) ** (J - 2 * a)
           * (u1 * w2 - w1 * u2) ** (J - 2 * b))
    pref = sqrt(Rational(fac(J + 1) * fac(J - 2 * a) * fac(J - 2 * b) * fac(J - 2 * c), 2 * c + 1))
    rhs = pref * _sum_abg(a, b, c, lambda al, be, ga:
        S.NegativeOne ** (a - b - ga)
        * u1 ** (a - al) * u2 ** (a + al) * v1 ** (b - be) * v2 ** (b + be) * w1 ** (c - ga) * w2 ** (c + ga)
        / sqrt(fac(a + al) * fac(a - al) * fac(b + be) * fac(b - be) * fac(c + ga) * fac(c - ga))
        * CG(a, b, c, al, be, -ga))
    return expand(lhs - rhs) == 0


def eq210c(a, b, c):
    # match the coefficient of x1^{J-2a} x2^{J-2b} x3^{J-2c} on both sides
    # (a polynomial identity in t1,t2,t3)
    J = a + b + c
    Delta = sqrt(Rational(fac(a + b - c) * fac(a - b + c) * fac(-a + b + c), fac(a + b + c + 1)))
    lhs = ((t2 - t3) ** (J - 2 * a) * (t3 - t1) ** (J - 2 * b) * (t1 - t2) ** (J - 2 * c)
           / (fac(J - 2 * a) * fac(J - 2 * b) * fac(J - 2 * c)))
    rhs = _sum_abg(a, b, c, lambda al, be, ga:
        S.NegativeOne ** (a - b - ga) / Delta
        * t1 ** (a + al) * t2 ** (b + be) * t3 ** (c + ga)
        / sqrt((2 * c + 1) * fac(a + al) * fac(a - al) * fac(b + be) * fac(b - be) * fac(c + ga) * fac(c - ga))
        * CG(a, b, c, al, be, -ga))
    return expand(lhs - rhs) == 0


def eq211(a, b, c, ga):
    J = a + b + c
    if a - b + ga < 0:
        return None
    pref = (sqrt(Rational(fac(J - 2 * b) * fac(c + ga) * (2 * c + 1),
                          fac(J - 2 * a) * fac(J - 2 * c) * fac(c - ga) * fac(J + 1)))
            / fac(a - b + ga))
    lhs = hyperexpand(pref * (t - 1) ** (J - 2 * c) * hyper((ga - c, a - b - c), (a - b + ga + 1,), t))
    rhs = S.Zero
    for be in rng(b):
        al = ga - be
        if abs(al) > a:
            continue
        rhs += (CG(a, b, c, al, be, ga) * t ** (b - be)
                / sqrt(fac(a + al) * fac(a - al) * fac(b + be) * fac(b - be)))
    return all(abs((lhs - rhs).subs(t, Rational(k, 7)).evalf()) < 1e-9 for k in (1, 3, 5, 9))


def eq212(a, b, c, ga):
    J = a + b + c
    gp = a - b                                      # gamma' = a - b
    pref = sqrt(Rational(2 * c + 1, fac(J + 1) * fac(J - 2 * c)))
    lhs = pref * Rotation.d(c, ga, gp, th).doit()
    rhs = S.Zero
    for be in rng(b):
        al = ga - be
        if abs(al) > a:
            continue
        rhs += (S.NegativeOne ** (b + be) * CG(a, b, c, al, be, ga)
                * cos(th / 2) ** (a + b + al - be) * sin(th / 2) ** (a + b - al + be)
                / sqrt(fac(a + al) * fac(a - al) * fac(b + be) * fac(b - be)))
    return all(abs((lhs - rhs).subs(th, Rational(k, 10)).evalf()) < 1e-9 for k in (3, 7, 11, 15))


def eq213(a, b, c):
    # match the coefficient of u^{J-2a} v^{J-2b} w^{J-2c} in 1/(1+u^2+v^2+w^2)
    J = a + b + c
    if J % 2:
        return None
    lco = (S.NegativeOne ** (J // 2) * fac(J // 2)
           / (fac((J - 2 * a) // 2) * fac((J - 2 * b) // 2) * fac((J - 2 * c) // 2)))
    rco = (S.NegativeOne ** (a - b) * sqrt(Rational(fac(J + 1), 2 * c + 1)) * CG(a, b, c, 0, 0, 0)
           / sqrt(fac(J - 2 * a) * fac(J - 2 * b) * fac(J - 2 * c)))
    return bool((lco - rco).simplify() == 0)


CASES = {
    "eq 8.8.2  binomials (1 var/mom)": (eq210, [(1, 1, 1), (1, 1, 2), (2, 1, 1),
                                              (Rational(1, 2), Rational(1, 2), 1)]),
    "eq 8.8.3 binomials (2 var/mom)": (eq210b, [(1, 1, 1), (1, 1, 2), (2, 1, 1)]),
    "eq 8.8.4 exponential":           (eq210c, [(1, 1, 1), (1, 1, 2), (2, 2, 2),
                                               (Rational(3, 2), 1, Rational(1, 2))]),
    "eq 8.8.5  hypergeometric":        (eq211, [(2, 1, 2, 1), (2, 2, 2, 0), (1, 1, 2, 0), (2, 1, 1, 0)]),
    "eq 8.8.6  Wigner d-function":     (eq212, [(1, 1, 1, 0), (1, 1, 2, 1), (2, 1, 2, -1),
                                              (Rational(1, 2), Rational(1, 2), 1, 0)]),
    "eq 8.8.7  Schwinger":             (eq213, [(1, 1, 2), (2, 2, 2), (1, 1, 0), (2, 1, 1), (2, 2, 0)]),
}


def run():
    print("Section 8.8 generating-function checks\n")
    all_ok = True
    for label, (fn, args) in CASES.items():
        results = [fn(*x) for x in args]
        results = [r for r in results if r is not None]
        ok = all(results)
        all_ok &= ok
        print(f"  [{'OK  ' if ok else 'FAIL'}] {label:32s} {sum(results)}/{len(results)} cases")
    print("\n  (not checked: eq 8.8.1, Regge determinant -- sum over all R-symbols, cumbersome)")
    print("\nALL CHECKED IDENTITIES HOLD" if all_ok else "SOME IDENTITIES FAILED")
    return all_ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
