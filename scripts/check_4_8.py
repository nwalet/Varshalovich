#!/usr/bin/env python3
r"""
Checks for Section 4.8 (Recursion relations for D^J) of Chapter 4, VMK.
Part 1: 4.8.1-4.8.9 (D^J <-> D^{J+-1}) and 4.8.10-4.8.15 (D^J <-> D^{J+-1/2}).

Each relation is  prefactor(a,b,g) * D^J_{M+dM, M'+dM'} = sum_k c_k D^{lvl_k}_{MM'}.
Coefficients are transcribed as printed; a FAIL flags an OCR/source error,
which is then corrected and re-checked.  D from wigner_d (0 if out of range).

Usage:  python3 check_4_8.py
"""
import os, sys, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy import Rational, S
from wigner_d import wigner_d, beta as _B

H = Rational(1, 2)
ANG = [(0.7, 1.1, 0.4), (1.3, 0.6, 2.0), (0.5, 2.1, 1.2)]
TOL = 1e-10
def sq(x):
    return cmath.sqrt(complex(x))

_dc = {}
def valid(J, M):
    return abs(M) <= J and (J - M) == int(J - M) and (J - M) >= 0

def dnum(J, M, N, b):
    key = (J, M, N)
    if key not in _dc:
        _dc[key] = wigner_d(J, M, N)
    return complex(_dc[key].subs(_B, float(b)).evalf(30))

def D(J, M, N, a, b, g):
    if J < 0 or not valid(J, M) or not valid(J, N):
        return 0j
    return cmath.exp(-1j*float(M)*a)*dnum(J, M, N, b)*cmath.exp(-1j*float(N)*g)

def f(x):
    return float(x)

def report(tag, worst):
    ok = worst < TOL
    print(f"  {tag:56s} {'PASS' if ok else 'FAIL'}  worst={worst:.1e}")
    return ok


def check(tag, Jset, Mset, pre, dM, dMp, coeffs):
    """coeffs(J,M,Mp) -> list of (c, level) for RHS c*D^level_{MM'}."""
    worst = 0.0
    for J in Jset:
        for M in Mset(J):
            for Mp in Mset(J):
                for (a, b, g) in ANG:
                    lhs = pre(a, b, g) * D(J, M+dM, Mp+dMp, a, b, g)
                    rhs = sum(c * D(lvl, M, Mp, a, b, g) for (c, lvl) in coeffs(J, M, Mp))
                    worst = max(worst, abs(lhs-rhs))
    return report(tag, worst)


def main():
    print("Section 4.8 recursions, part 1\n")
    Jint = [S(1), Rational(3,2), S(2), Rational(5,2)]
    full = lambda J: [J - k for k in range(int(2*J)+1)]
    import math
    cosb = lambda a,b,g: math.cos(b)
    sinb = lambda a,b,g: math.sin(b)
    ok = True

    # 4.8.1  cos b D^J = ...
    ok &= check("4.8.1  cos b D^J", Jint, full, cosb, 0, 0, lambda J,M,Mp: [
        (sq((J**2-M**2)*(J**2-Mp**2))/(J*(2*J+1)), J-1),
        (M*Mp/(J*(J+1)), J),
        (sq(((J+1)**2-M**2)*((J+1)**2-Mp**2))/((J+1)*(2*J+1)), J+1)])
    # 4.8.2  sin b e^{ia} D^J_{M+1,M'}
    ok &= check("4.8.2  sin b e^{ia} D^J_{M+1,M'}", Jint, full,
                lambda a,b,g: math.sin(b)*cmath.exp(1j*a), 1, 0, lambda J,M,Mp: [
        (-sq((J+M)*(J+M+1)*(J**2-Mp**2))/(J*(2*J+1)), J-1),
        (Mp*sq((J-M)*(J+M+1))/(J*(J+1)), J),
        (sq((J-M)*(J-M+1)*((J+1)**2-Mp**2))/((J+1)*(2*J+1)), J+1)])
    # 4.8.3  sin b e^{-ia} D^J_{M-1,M'}
    ok &= check("4.8.3  sin b e^{-ia} D^J_{M-1,M'}", Jint, full,
                lambda a,b,g: math.sin(b)*cmath.exp(-1j*a), -1, 0, lambda J,M,Mp: [
        (sq((J-M)*(J-M+1)*(J**2-Mp**2))/(J*(2*J+1)), J-1),
        (Mp*sq((J+M)*(J-M+1))/(J*(J+1)), J),
        (-sq((J+M)*(J+M+1)*((J+1)**2-Mp**2))/((J+1)*(2*J+1)), J+1)])
    # 4.8.4  sin b e^{ig} D^J_{M,M'+1}
    ok &= check("4.8.4  sin b e^{ig} D^J_{M,M'+1}", Jint, full,
                lambda a,b,g: math.sin(b)*cmath.exp(1j*g), 0, 1, lambda J,M,Mp: [
        (sq((J**2-M**2)*(J+Mp)*(J+Mp+1))/(J*(2*J+1)), J-1),
        (-M*sq((J-Mp)*(J+Mp+1))/(J*(J+1)), J),
        (-sq(((J+1)**2-M**2)*(J-Mp)*(J-Mp+1))/((J+1)*(2*J+1)), J+1)])
    # 4.8.5  sin b e^{-ig} D^J_{M,M'-1}
    ok &= check("4.8.5  sin b e^{-ig} D^J_{M,M'-1}", Jint, full,
                lambda a,b,g: math.sin(b)*cmath.exp(-1j*g), 0, -1, lambda J,M,Mp: [
        (-sq((J**2-M**2)*(J-Mp)*(J-Mp+1))/(J*(2*J+1)), J-1),
        (-M*sq((J+Mp)*(J-Mp+1))/(J*(J+1)), J),
        (sq(((J+1)**2-M**2)*(J+Mp)*(J+Mp+1))/((J+1)*(2*J+1)), J+1)])
    # 4.8.6  (1+cos b) e^{i(a+g)} D^J_{M+1,M'+1}
    ok &= check("4.8.6  (1+cosb)e^{i(a+g)} D^J_{M+1,M'+1}", Jint, full,
                lambda a,b,g: (1+math.cos(b))*cmath.exp(1j*(a+g)), 1, 1, lambda J,M,Mp: [
        (sq((J+M+1)*(J+M)*(J+Mp+1)*(J+Mp))/(J*(2*J+1)), J-1),
        (sq((J-M)*(J+M+1)*(J-Mp)*(J+Mp+1))/(J*(J+1)), J),
        (sq((J-M)*(J-M+1)*(J-Mp)*(J-Mp+1))/((J+1)*(2*J+1)), J+1)])
    # 4.8.7  (1+cos b) e^{-i(a+g)} D^J_{M-1,M'-1}   -- last coeff corrected
    ok &= check("4.8.7  (1+cosb)e^{-i(a+g)} [c3 corrected]", Jint, full,
                lambda a,b,g: (1+math.cos(b))*cmath.exp(-1j*(a+g)), -1, -1, lambda J,M,Mp: [
        (sq((J-M)*(J-M+1)*(J-Mp)*(J-Mp+1))/(J*(2*J+1)), J-1),
        (sq((J+M)*(J-M+1)*(J+Mp)*(J-Mp+1))/(J*(J+1)), J),
        (sq((J+M)*(J+M+1)*(J+Mp)*(J+Mp+1))/((J+1)*(2*J+1)), J+1)])
    # 4.8.8  (1-cos b) e^{i(a-g)} D^J_{M+1,M'-1}
    ok &= check("4.8.8  (1-cosb)e^{i(a-g)} D^J_{M+1,M'-1}", Jint, full,
                lambda a,b,g: (1-math.cos(b))*cmath.exp(1j*(a-g)), 1, -1, lambda J,M,Mp: [
        (sq((J+M)*(J+M+1)*(J-Mp)*(J-Mp+1))/(J*(2*J+1)), J-1),
        (-sq((J-M)*(J+M+1)*(J+Mp)*(J-Mp+1))/(J*(J+1)), J),
        (sq((J-M)*(J-M+1)*(J+Mp)*(J+Mp+1))/((J+1)*(2*J+1)), J+1)])
    # 4.8.9  (1-cos b) e^{-i(a-g)} D^J_{M-1,M'+1}
    ok &= check("4.8.9  (1-cosb)e^{-i(a-g)} D^J_{M-1,M'+1}", Jint, full,
                lambda a,b,g: (1-math.cos(b))*cmath.exp(-1j*(a-g)), -1, 1, lambda J,M,Mp: [
        (sq((J-M)*(J-M+1)*(J+Mp)*(J+Mp+1))/(J*(2*J+1)), J-1),
        (-sq((J+M)*(J-M+1)*(J-Mp)*(J+Mp+1))/(J*(J+1)), J),
        (sq((J+M)*(J+M+1)*(J-Mp)*(J-Mp+1))/((J+1)*(2*J+1)), J+1)])

    # ---- 4.8.10-4.8.13: D^{J+-1/2}.  M ranges over projections of level J-1/2 ----
    half = lambda J: [J - H - k for k in range(int(2*(J-H))+1)]
    print()
    # phases as PRINTED are e^{i(a+a)/2} etc.; test the corrected e^{i(a+g)/2}
    ok &= check("4.8.10 cos(b/2)e^{i(a+g)/2} D^J_{M+1/2,M'+1/2}", Jint, half,
                lambda a,b,g: math.cos(b/2)*cmath.exp(1j*(a+g)/2), H, H, lambda J,M,Mp: [
        (sq((J+M+H)*(J+Mp+H))/(2*J+1), J-H),
        (sq((J-M+H)*(J-Mp+H))/(2*J+1), J+H)])
    ok &= check("4.8.11 sin(b/2)e^{i(a-g)/2} D^J_{M+1/2,M'-1/2}", Jint, half,
                lambda a,b,g: math.sin(b/2)*cmath.exp(1j*(a-g)/2), H, -H, lambda J,M,Mp: [
        (-sq((J+M+H)*(J-Mp+H))/(2*J+1), J-H),
        (sq((J-M+H)*(J+Mp+H))/(2*J+1), J+H)])
    ok &= check("4.8.12 cos(b/2)e^{-i(a+g)/2} D^J_{M-1/2,M'-1/2}", Jint, half,
                lambda a,b,g: math.cos(b/2)*cmath.exp(-1j*(a+g)/2), -H, -H, lambda J,M,Mp: [
        (sq((J-M+H)*(J-Mp+H))/(2*J+1), J-H),
        (sq((J+M+H)*(J+Mp+H))/(2*J+1), J+H)])
    ok &= check("4.8.13 sin(b/2)e^{-i(a-g)/2} D^J_{M-1/2,M'+1/2}", Jint, half,
                lambda a,b,g: math.sin(b/2)*cmath.exp(-1j*(a-g)/2), -H, H, lambda J,M,Mp: [
        (sq((J-M+H)*(J+Mp+H))/(2*J+1), J-H),
        (-sq((J+M+H)*(J-Mp+H))/(2*J+1), J+H)])

    # ---- 4.8.14, 4.8.15: D^J in terms of D^{J-1/2} (phases x->gamma corrected) ----
    print()
    def check_inv(tag, cond, terms):
        worst = 0.0
        for J in Jint:
            for M in full(J):
                for Mp in full(J):
                    if cond(J, M, Mp):
                        continue
                    for (a, b, g) in ANG:
                        lhs = D(J, M, Mp, a, b, g)
                        rhs = sum(t(J, M, Mp, a, b, g) for t in terms)
                        worst = max(worst, abs(lhs-rhs))
        return report(tag, worst)
    # 4.8.14 (M'!=J)
    ok &= check_inv("4.8.14 D^J = ... D^{J-1/2}  (M'!=J)", lambda J,M,Mp: Mp == J, [
        lambda J,M,Mp,a,b,g: sq((J-M)/(J-Mp))*math.cos(b/2)*cmath.exp(1j*(a+g)/2)*D(J-H,M+H,Mp+H,a,b,g),
        lambda J,M,Mp,a,b,g: -sq((J+M)/(J-Mp))*math.sin(b/2)*cmath.exp(-1j*(a-g)/2)*D(J-H,M-H,Mp+H,a,b,g)])
    # 4.8.15 (M'!=-J)
    ok &= check_inv("4.8.15 D^J = ... D^{J-1/2}  (M'!=-J)", lambda J,M,Mp: Mp == -J, [
        lambda J,M,Mp,a,b,g: sq((J-M)/(J+Mp))*math.sin(b/2)*cmath.exp(1j*(a-g)/2)*D(J-H,M+H,Mp-H,a,b,g),
        lambda J,M,Mp,a,b,g: sq((J+M)/(J+Mp))*math.cos(b/2)*cmath.exp(-1j*(a+g)/2)*D(J-H,M-H,Mp-H,a,b,g)])

    # ---- 4.8.16-4.8.21: relations among D^J at the SAME level (shifted M/M') ----
    print()
    def check_same(tag, pre, dM, dMp, terms):
        worst = 0.0
        for J in Jint:
            for M in full(J):
                for Mp in full(J):
                    for (a, b, g) in ANG:
                        try:
                            lhs = pre(J,M,Mp,a,b,g) * D(J, M+dM, Mp+dMp, a, b, g)
                            rhs = sum(t(J,M,Mp,a,b,g) for t in terms)
                            d = abs(complex(lhs) - complex(rhs))
                        except (ZeroDivisionError, TypeError, ValueError):
                            continue
                        if math.isfinite(d):
                            worst = max(worst, d)
        return report(tag, worst)
    e = cmath.exp
    ok &= check_same("4.8.16 (-M+M'cosb)/sinb D^J", lambda J,M,Mp,a,b,g:(-M+Mp*math.cos(b))/math.sin(b), 0,0, [
        lambda J,M,Mp,a,b,g: 0.5*sq((J+Mp)*(J-Mp+1))*e(-1j*g)*D(J,M,Mp-1,a,b,g),
        lambda J,M,Mp,a,b,g: 0.5*sq((J-Mp)*(J+Mp+1))*e( 1j*g)*D(J,M,Mp+1,a,b,g)])
    ok &= check_same("4.8.17 (M'-Mcosb)/sinb D^J", lambda J,M,Mp,a,b,g:(Mp-M*math.cos(b))/math.sin(b), 0,0, [
        lambda J,M,Mp,a,b,g: 0.5*sq((J+M)*(J-M+1))*e(-1j*a)*D(J,M-1,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: 0.5*sq((J-M)*(J+M+1))*e( 1j*a)*D(J,M+1,Mp,a,b,g)])
    # 4.8.18 (from scan): B denom (J-M)(J+M+1); C denom (J-M)(J+M+1)
    ok &= check_same("4.8.18 D^J_{M+1,M'} e^{ia}", lambda J,M,Mp,a,b,g:e(1j*a), 1,0, [
        lambda J,M,Mp,a,b,g: sq((J+Mp)*(J-Mp+1)/((J-M)*(J+M+1)))*(1+math.cos(b))/2*e(-1j*g)*D(J,M,Mp-1,a,b,g),
        lambda J,M,Mp,a,b,g: Mp*math.sin(b)/sq((J-M)*(J+M+1))*D(J,M,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: -sq((J-Mp)*(J+Mp+1)/((J-M)*(J+M+1)))*(1-math.cos(b))/2*e(1j*g)*D(J,M,Mp+1,a,b,g)])
    # 4.8.19 (from scan): all denoms (J+M)(J-M+1)
    ok &= check_same("4.8.19 D^J_{M-1,M'} e^{-ia}", lambda J,M,Mp,a,b,g:e(-1j*a), -1,0, [
        lambda J,M,Mp,a,b,g: -sq((J+Mp)*(J-Mp+1)/((J+M)*(J-M+1)))*(1-math.cos(b))/2*e(-1j*g)*D(J,M,Mp-1,a,b,g),
        lambda J,M,Mp,a,b,g: Mp*math.sin(b)/sq((J+M)*(J-M+1))*D(J,M,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: sq((J-Mp)*(J+Mp+1)/((J+M)*(J-M+1)))*(1+math.cos(b))/2*e(1j*g)*D(J,M,Mp+1,a,b,g)])
    # 4.8.20 (from scan): diagonal uses M (not M'); C num (J-M)(J+M+1)
    ok &= check_same("4.8.20 D^J_{M,M'+1} e^{ig}", lambda J,M,Mp,a,b,g:e(1j*g), 0,1, [
        lambda J,M,Mp,a,b,g: sq((J+M)*(J-M+1)/((J-Mp)*(J+Mp+1)))*(1+math.cos(b))/2*e(-1j*a)*D(J,M-1,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: -M*math.sin(b)/sq((J-Mp)*(J+Mp+1))*D(J,M,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: -sq((J-M)*(J+M+1)/((J-Mp)*(J+Mp+1)))*(1-math.cos(b))/2*e(1j*a)*D(J,M+1,Mp,a,b,g)])
    # 4.8.21 (from scan): diagonal uses M (not M'); C num (J-M)(J+M+1)
    ok &= check_same("4.8.21 D^J_{M,M'-1} e^{-ig}", lambda J,M,Mp,a,b,g:e(-1j*g), 0,-1, [
        lambda J,M,Mp,a,b,g: -sq((J+M)*(J-M+1)/((J+Mp)*(J-Mp+1)))*(1-math.cos(b))/2*e(-1j*a)*D(J,M-1,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: -M*math.sin(b)/sq((J+Mp)*(J-Mp+1))*D(J,M,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: sq((J-M)*(J+M+1)/((J+Mp)*(J-Mp+1)))*(1+math.cos(b))/2*e(1j*a)*D(J,M+1,Mp,a,b,g)])

    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
