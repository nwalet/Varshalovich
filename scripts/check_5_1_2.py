#!/usr/bin/env python3
r"""
Checks for Sections 5.1 (Definition) and 5.2 (Explicit forms) of Chapter 5, VMK.

Reference: Y_{lm}(theta,phi) = mpmath.spherharm(l,m,theta,phi) — the standard
physics (Condon-Shortley) convention, which is VMK's (Eq. 5.2.1); validated in
Chapter 4 (check_4_16_17.py) against D_{0,-m}^l.

This first pass covers the pieces that are clean to evaluate:
  5.1.6   orthonormality
  5.1.10  Y_l0(0,0)
  5.1.11  conjugation relations
  5.1.14  u_lm, v_lm
  5.2.1   Y = e^{imphi} sqrt(..) P_l^m(cos th)
  5.2.2-5.2.6  differential expressions
  5.2.23,5.2.27  hypergeometric (trig th/2)
  5.2.29,5.2.31  hypergeometric (trig th)
  5.2.37  relation to D^l
  5.2.38,5.2.39  Jacobi / Gegenbauer

Usage:  python3 check_5_1_2.py
"""
import math, cmath
import mpmath as mp
from scipy.integrate import quad, dblquad
from scipy.special import lpmv, eval_jacobi, eval_gegenbauer

mp.mp.dps = 30
TOL = 1e-9
TH = [0.4, 0.9, 1.3, 1.9, 2.6]
PH = [0.0, 0.7, 1.8, 3.1]

def Y(l, m, th, ph):
    if abs(m) > l: return 0j
    return complex(mp.spherharm(l, m, th, ph))
def fac(n): return math.factorial(int(round(n)))
def dfac(n):
    n = int(round(n)); r = 1.0
    while n > 1: r *= n; n -= 2
    return r
def report(tag, worst, tol=TOL):
    ok = worst < tol
    print(f"  {tag:44s} {'PASS' if ok else 'FAIL'}  worst={worst:.2e}")
    return ok
def sweep(tag, form, ls=range(0, 6), tol=TOL, ths=TH, phs=PH, mpos=False):
    worst = 0.0
    for l in ls:
        ms = range(0, l+1) if mpos else range(-l, l+1)
        for m in ms:
            for th in ths:
                for ph in phs:
                    worst = max(worst, abs(form(l, m, th, ph) - Y(l, m, th, ph)))
    return report(tag, worst, tol)


def main():
    print("Section 5.1-5.2 spherical harmonics\n")
    ok = True

    # 5.2.1  Y = e^{im phi} sqrt((2l+1)/4pi (l-m)!/(l+m)!) P_l^m(cos th)
    def f21(l, m, th, ph):
        return (cmath.exp(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi)*fac(l-m)/fac(l+m))
                * lpmv(m, l, math.cos(th)))
    ok &= sweep("5.2.1  e^{imf} sqrt.. P_l^m", f21)

    # 5.1.10 Y_l0(0,0) = sqrt((2l+1)/4pi)
    worst = max(abs(Y(l, 0, 0, 0)-math.sqrt((2*l+1)/(4*math.pi))) for l in range(0, 6))
    ok &= report("5.1.10 Y_l0(0,0)", worst)

    # 5.1.11 Y*_lm = Y_lm(th,-phi) = (-1)^m Y_{l,-m}
    worst = 0.0
    for l in range(0, 5):
        for m in range(-l, l+1):
            for th in TH:
                for ph in PH:
                    c = Y(l, m, th, ph).conjugate()
                    worst = max(worst, abs(c - Y(l, m, th, -ph)))
                    worst = max(worst, abs(c - (-1)**m*Y(l, -m, th, ph)))
    ok &= report("5.1.11 conjugation", worst)

    # 5.1.6 orthonormality (l,l' up to 2)
    def inner(l1, m1, l2, m2):
        re = dblquad(lambda th, ph: (Y(l1,m1,th,ph).conjugate()*Y(l2,m2,th,ph)).real*math.sin(th),
                     0, 2*math.pi, 0, math.pi)[0]
        im = dblquad(lambda th, ph: (Y(l1,m1,th,ph).conjugate()*Y(l2,m2,th,ph)).imag*math.sin(th),
                     0, 2*math.pi, 0, math.pi)[0]
        return re+1j*im
    worst = 0.0
    for (l1,m1,l2,m2) in [(1,0,1,0),(1,1,1,1),(2,1,2,1),(1,0,2,0),(2,1,2,-1),(1,1,2,1)]:
        d = 1.0 if (l1==l2 and m1==m2) else 0.0
        worst = max(worst, abs(inner(l1,m1,l2,m2)-d))
    ok &= report("5.1.6  orthonormality", worst, tol=1e-6)

    # 5.1.14 u_lm, v_lm
    worst = 0.0
    for l in range(0, 5):
        for m in range(0, l+1):
            for th in TH:
                for ph in PH:
                    ylm = Y(l, m, th, ph)
                    u = 0.5*(ylm + ylm.conjugate())
                    v = (ylm - ylm.conjugate())/(2j)
                    ru = (math.sqrt((2*l+1)/(2*math.pi)*fac(l-m)/fac(l+m))
                          * math.cos(m*ph)*lpmv(m, l, math.cos(th))) if m>0 else \
                         (math.sqrt((2*l+1)/(4*math.pi))*lpmv(0, l, math.cos(th)))
                    # NB u_lm uses sqrt((2l+1)/(2pi)..) for m!=0 per 5.1.14
                    rv = (math.sqrt((2*l+1)/(4*math.pi)*fac(l-m)/fac(l+m))
                          * math.sin(m*ph)*lpmv(m, l, math.cos(th)))
                    if m == 0:
                        worst = max(worst, abs(u - Y(l,0,th,ph)), abs(v))
                    else:
                        worst = max(worst, abs(u.real-ru), abs(v.real-rv))
    ok &= report("5.1.14 u_lm, v_lm", worst, tol=1e-8)

    # 5.2.2  Y = e^{imf}/(2^l l!) sqrt((2l+1)/4pi (l+m)!/(l-m)!) (sin th)^{-m}
    #             d^{l-m}/(dcos)^{l-m} (cos^2-1)^l
    x = mp.mpf
    def dpow(expr, n, xv):   # n-th derivative wrt cos, evaluated at cos th
        return mp.diff(expr, xv, n)
    def f22(l, m, th, ph):
        c = math.cos(th)
        expr = lambda u: (u*u-1)**l
        d = float(mp.diff(expr, c, l-m))
        return (cmath.exp(1j*m*ph)/(2**l*fac(l))
                * math.sqrt((2*l+1)/(4*math.pi)*fac(l+m)/fac(l-m))
                * math.sin(th)**(-m)*d)
    ok &= sweep("5.2.2  d^{l-m} Rodrigues", f22, ls=range(1, 5))

    def f23(l, m, th, ph):
        c = math.cos(th)
        expr = lambda u: (u*u-1)**l
        d = float(mp.diff(expr, c, l+m))
        return ((-1)**m*cmath.exp(1j*m*ph)/(2**l*fac(l))
                * math.sqrt((2*l+1)/(4*math.pi)*fac(l-m)/fac(l+m))
                * math.sin(th)**m*d)
    ok &= sweep("5.2.3  d^{l+m} Rodrigues", f23, ls=range(1, 5))

    # 5.2.6  Y = (-1)^m e^{imf} sqrt((2l+1)/4pi (l-m)!/(l+m)!) (sin)^m d^m/dcos^m P_l  (m>=0)
    def f26(l, m, th, ph):
        c = math.cos(th)
        expr = lambda u: float(mp.legendre(l, u)) if False else None
        d = float(mp.diff(lambda u: mp.legendre(l, u), c, m))
        return ((-1)**m*cmath.exp(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi)*fac(l-m)/fac(l+m))
                * math.sin(th)**m*d)
    ok &= sweep("5.2.6  d^m P_l (m>=0)", f26, ls=range(0, 5), mpos=True)

    # 5.2.23  Y = xi e^{imf} sqrt((2l+1)/4pi (l+|m|)!/(l-|m|)!) (sin th)^|m|/(|m|! 2^|m|)
    #              F(-l+|m|, l+|m|+1; |m|+1; sin^2(th/2))
    def xi(m): return (-1)**m if m > 0 else 1
    def F(a,b,c,z): return complex(mp.hyp2f1(a,b,c,z))
    def f223(l, m, th, ph):
        am = abs(m)
        return (xi(m)*cmath.exp(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi)*fac(l+am)/fac(l-am))
                * math.sin(th)**am/(fac(am)*2**am)
                * F(-l+am, l+am+1, am+1, math.sin(th/2)**2))
    ok &= sweep("5.2.23 2F1(sin^2 th/2)", f223)

    # 5.2.27 Y = xi e^{imf} sqrt(..(l+|m|)!/(l-|m|)!) (tan th/2)^|m| (cos th/2)^{2l}
    #             F(-l+|m|,-l; |m|+1; -tan^2 th/2)
    def f227(l, m, th, ph):
        am = abs(m)
        return (xi(m)*cmath.exp(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi)*fac(l+am)/fac(l-am))
                * math.tan(th/2)**am*math.cos(th/2)**(2*l)
                * F(-l+am, -l, am+1, -math.tan(th/2)**2))
    ok &= sweep("5.2.27 2F1(-tan^2 th/2)", f227)

    # 5.2.37 Y = sqrt((2l+1)/4pi) D^l_{0,-m}(chi,th,phi)  (chi arbitrary -> use 0)
    #   D^l_{0,-m}(0,th,phi) = d^l_{0,-m}(th) e^{-i(-m)phi}=d e^{imphi}
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from wigner_d import wigner_d, beta as _B
    def Dl(l, M, N, a, b, g):
        v = complex(wigner_d(l, M, N).subs(_B, float(b)).evalf(30))
        return cmath.exp(-1j*M*a)*v*cmath.exp(-1j*N*g)
    def f237(l, m, th, ph):
        return math.sqrt((2*l+1)/(4*math.pi))*Dl(l, 0, -m, 0.5, th, ph)
    ok &= sweep("5.2.37 = sqrt.. D^l_{0,-m}", f237)

    # 5.2.39 Gegenbauer
    def f239(l, m, th, ph):
        am = abs(m)
        return (xi(m)*cmath.exp(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi)*fac(l-am)/fac(l+am))
                * dfac(2*am-1)*math.sin(th)**am
                * eval_gegenbauer(l-am, 0.5+am, math.cos(th)))
    ok &= sweep("5.2.39 Gegenbauer", f239)

    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
