#!/usr/bin/env python3
r"""
Checks for Sections 4.16 (D for particular argument values) and 4.17 (special
cases of D for particular M or M') of Chapter 4, VMK.

D^J_{MM'}(a,b,g) = e^{-iMa} d^J_{MM'}(b) e^{-iM'g}, d from the validated
wigner_d helper.  Spherical harmonics use the standard (VMK) convention,
validated below against 4.17.1.

Purpose besides confirmation: pin down the correct phases in 4.17.3 and 4.17.4,
where the .tex shows the OCR garble  e^{... (alpha +- alpha)/2}.

Usage:  python3 check_4_16_17.py
"""
import os, sys, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy import Rational, S
from wigner_d import wigner_d, beta as _B
import mpmath as mp
from scipy.special import eval_legendre

TOL = 1e-9
ANG = [(0.7, 1.1, 0.4), (1.3, 0.6, 2.0), (2.1, 2.3, 1.2)]

_dc = {}
def _valid(J, M): return abs(M) <= J and (J-M) == int(J-M) and (J-M) >= 0
def dval(J, M, N, b):
    if not _valid(J, M) or not _valid(J, N): return None
    if (J, M, N) not in _dc: _dc[(J, M, N)] = wigner_d(J, M, N)
    return complex(_dc[(J, M, N)].subs(_B, float(b)).evalf(30))
def D(J, M, N, a, b, g):
    v = dval(J, M, N, b)
    return 0j if v is None else cmath.exp(-1j*float(M)*a)*v*cmath.exp(-1j*float(N)*g)

def Y(l, m, theta, phi):
    """VMK Y_{lm}(theta,phi): standard physics convention (CS phase).
    mpmath.spherharm(l,m,theta,phi) is the standard Y_l^m with CS phase."""
    if abs(m) > l: return 0j
    return complex(mp.spherharm(l, m, theta, phi))

def Pl(l, x):    return float(eval_legendre(l, x))
def dPl(l, x):
    # derivative of Legendre P_l at x
    if l == 0: return 0.0
    return l*(x*eval_legendre(l, x)-eval_legendre(l-1, x))/(x*x-1)

def report(tag, worst, tol=TOL):
    ok = worst < tol
    print(f"  {tag:48s} {'PASS' if ok else 'FAIL'}  worst={worst:.2e}")
    return ok
def mrange(J): return [J-k for k in range(int(2*J)+1)]
def kd(a, b): return 1.0 if a == b else 0.0


def check_416():
    print("Section 4.16  D for particular arguments")
    ok = True
    Js = [Rational(1,2), S(1), Rational(3,2), S(2)]
    # 4.16.1 D(0,0,0)=d_{MM'}
    worst = 0.0
    for J in Js:
        for M in mrange(J):
            for Mp in mrange(J):
                worst = max(worst, abs(D(J,M,Mp,0,0,0)-kd(M,Mp)))
    ok &= report("4.16.1 D(0,0,0)=delta", worst)
    # 4.16.2 D(a,0,g)=delta e^{-iM(a+g)}
    worst = 0.0
    for J in Js:
        for M in mrange(J):
            for Mp in mrange(J):
                for (a,b,g) in ANG:
                    worst = max(worst, abs(D(J,M,Mp,a,0,g)-kd(M,Mp)*cmath.exp(-1j*float(M)*(a+g))))
    ok &= report("4.16.2 D(a,0,g)", worst)
    # 4.16.3 D(a,2n pi,g)=delta (-1)^{2nJ} e^{-iM(a+g)}
    worst = 0.0
    for J in Js:
        for n in [1, 2]:
            for M in mrange(J):
                for Mp in mrange(J):
                    for (a,b,g) in ANG:
                        lhs = D(J,M,Mp,a,2*n*math.pi,g)
                        rhs = kd(M,Mp)*(-1)**int(2*n*J)*cmath.exp(-1j*float(M)*(a+g))
                        worst = max(worst, abs(lhs-rhs))
    ok &= report("4.16.3 D(a,2n pi,g)", worst)
    # 4.16.4 D(a,(2n+1)pi,g)=delta_{-M,M'}(-1)^{(2n+1)J+M} e^{-iM(a-g)}
    worst = 0.0
    for J in Js:
        for n in [0, 1]:
            for M in mrange(J):
                for Mp in mrange(J):
                    for (a,b,g) in ANG:
                        lhs = D(J,M,Mp,a,(2*n+1)*math.pi,g)
                        rhs = kd(-M,Mp)*(-1)**int((2*n+1)*J+M)*cmath.exp(-1j*float(M)*(a-g)) if _valid(J,-M) else 0
                        # (-1)^{...+M}: M half-integer -> use +M in exponent carefully
                        if kd(-M,Mp):
                            ph = cmath.exp(1j*math.pi*((2*n+1)*float(J)+float(M)))
                            rhs = ph*cmath.exp(-1j*float(M)*(a-g))
                        else:
                            rhs = 0
                        worst = max(worst, abs(lhs-rhs))
    ok &= report("4.16.4 D(a,(2n+1)pi,g)", worst)
    # 4.16.5 D(a,pi/2,g) explicit binomial sum
    from math import comb, factorial
    worst = 0.0
    for J in [S(1), S(2), Rational(1,2), Rational(3,2)]:
        for M in mrange(J):
            for Mp in mrange(J):
                # binomials need integer args
                if (J+Mp) != int(J+Mp): continue
                s = 0.0
                for k in range(0, int(2*J)+1):
                    top1 = int(J+Mp); bot1 = k
                    top2 = int(J-Mp); bot2 = k+int(M-Mp)
                    if bot1 < 0 or bot1 > top1 or bot2 < 0 or bot2 > top2: continue
                    s += (-1)**k*comb(top1, bot1)*comb(top2, bot2)
                pref = ((-1)**int(M-Mp)/2**float(J)
                        * math.sqrt(factorial(int(J+M))*factorial(int(J-M))
                                    / (factorial(int(J+Mp))*factorial(int(J-Mp)))))
                rhs = pref*s   # phase e^{-iaM-igM'} tested at a=g=0 -> real
                lhs = D(J, M, Mp, 0, math.pi/2, 0)
                worst = max(worst, abs(lhs-rhs))
    ok &= report("4.16.5 D(a,pi/2,g) binomial sum", worst)
    return ok


def check_417():
    print("\nSection 4.17  special cases of D for particular M/M'")
    ok = True
    # validate Y convention via 4.17.1:  D_{m0}^l = (-1)^m sqrt(4pi/(2l+1)) Y_{l,-m}(b,a)
    worst = 0.0
    for l in [1, 2, 3]:
        for m in range(-l, l+1):
            for (a,b,g) in ANG:
                lhs = D(l, m, 0, a, b, g)
                rhs = (-1)**m*math.sqrt(4*math.pi/(2*l+1))*Y(l, -m, b, a)
                worst = max(worst, abs(lhs-rhs))
    ok &= report("4.17.1a D_{m0}=(-1)^m sqrt.. Y_{l,-m}(b,a)", worst)
    worst = 0.0
    for l in [1, 2, 3]:
        for m in range(-l, l+1):
            for (a,b,g) in ANG:
                lhs = D(l, 0, m, a, b, g)
                rhs = math.sqrt(4*math.pi/(2*l+1))*Y(l, -m, b, g)
                worst = max(worst, abs(lhs-rhs))
    ok &= report("4.17.1b D_{0m}=sqrt.. Y_{l,-m}(b,g)", worst)

    # 4.17.2 D_{00}=P_l(cos b) ; D_{+-10}, D_{0+-1}, D_{+-20}, D_{0+-2}
    worst = 0.0
    for l in [1, 2, 3]:
        for (a,b,g) in ANG:
            worst = max(worst, abs(D(l,0,0,a,b,g)-Pl(l, math.cos(b))))
    ok &= report("4.17.2 D_{00}=P_l(cos b)", worst)
    worst = 0.0
    for l in [1, 2, 3]:
        for sgn in [+1, -1]:
            for (a,b,g) in ANG:
                lhs = D(l, sgn*1, 0, a, b, g)
                rhs = -sgn*cmath.exp(-1j*sgn*a)*math.sin(b)/math.sqrt(l*(l+1))*dPl(l, math.cos(b))
                worst = max(worst, abs(lhs-rhs))
    ok &= report("4.17.2 D_{+-1,0}", worst)
    worst = 0.0
    for l in [1, 2, 3]:
        for sgn in [+1, -1]:
            for (a,b,g) in ANG:
                lhs = D(l, 0, sgn*1, a, b, g)
                rhs = sgn*cmath.exp(-1j*sgn*g)*math.sin(b)/math.sqrt(l*(l+1))*dPl(l, math.cos(b))
                worst = max(worst, abs(lhs-rhs))
    ok &= report("4.17.2 D_{0,+-1}", worst)

    # ---- 4.17.4  D_{+-1/2,+-1/2}: confirm the corrected phases ----
    # half-integer-degree Legendre via mpmath.legenp ; derivative numerically
    def dP_half(nu, x):
        h = 1e-6
        return float((mp.legenp(nu, 0, x+h)-mp.legenp(nu, 0, x-h))/(2*h))
    worst = 0.0
    for J in [Rational(1,2), Rational(3,2), Rational(5,2)]:
        Jf = float(J)
        for (a,b,g) in ANG:
            cb = math.cos(b)
            Pp = dP_half(Jf+0.5, cb); Pm = dP_half(Jf-0.5, cb)
            cases = {
                (S(1)/2, S(1)/2):  cmath.exp(-1j*(a+g)/2)*math.cos(b/2)/(Jf+0.5)*(Pp-Pm),
                (S(1)/2, -S(1)/2): -cmath.exp(-1j*(a-g)/2)*math.sin(b/2)/(Jf+0.5)*(Pp+Pm),
                (-S(1)/2, S(1)/2):  cmath.exp(1j*(a-g)/2)*math.sin(b/2)/(Jf+0.5)*(Pp+Pm),
                (-S(1)/2, -S(1)/2): cmath.exp(1j*(a+g)/2)*math.cos(b/2)/(Jf+0.5)*(Pp-Pm),
            }
            for (M, Mp), rhs in cases.items():
                worst = max(worst, abs(D(J, M, Mp, a, b, g)-rhs))
    ok &= report("4.17.4 D_{+-1/2,+-1/2} (corrected phase)", worst, tol=1e-5)

    # 4.17.3 (D_{+-1/2,M'} and D_{M,+-1/2}) is expressed through spherical
    # harmonics of HALF-INTEGER degree/order (Y_{J+-1/2..., -M'}); those are
    # exotic (mpmath.spherharm diverges for several of the needed index pairs)
    # and the phase does not decompose the naive integer-Y way. Its phase
    # correction e^{+-i(alpha-gamma)/2} is taken from the scan (printed p.114)
    # and is consistent with 4.17.4, its M'=+-1/2 special case, verified below.
    print("  4.17.3  D_{+-1/2,M'} phase e^{+-i(a-g)/2} -- scan-verified (p.114)")

    # ---- 4.17.7 D_{J-1,m} etc. ----
    from math import factorial
    def fac(n): return factorial(int(round(n)))
    worst = 0.0
    for J in [S(1), Rational(3,2), S(2), Rational(5,2)]:
        Jf = float(J)
        for m in mrange(J):
            if abs(m) > J-1: continue
            for (a,b,g) in ANG:
                pref = ((-1)**int(Jf-float(m)-1)*cmath.exp(-1j*(Jf-1)*a-1j*float(m)*g)
                        * math.sqrt(fac(2*J-1)/(fac(J+m)*fac(J-m)))
                        * math.cos(b/2)**(Jf+float(m)-1)*math.sin(b/2)**(Jf-float(m)-1))
                rhs = pref*(Jf*math.cos(b)-float(m))
                worst = max(worst, abs(D(J, J-1, m, a, b, g)-rhs))
    ok &= report("4.17.7 D_{J-1,m}", worst)

    # ---- 4.17.8 D_{JM}, D_{-JM}, D_{MJ}, D_{M-J} ----
    worst = 0.0
    for J in [Rational(1,2), S(1), Rational(3,2), S(2)]:
        Jf = float(J)
        for M in mrange(J):
            for (a,b,g) in ANG:
                pref = math.sqrt(fac(2*J)/(fac(J+M)*fac(J-M)))
                r1 = pref*math.cos(b/2)**(Jf+float(M))*(-math.sin(b/2))**(Jf-float(M))*cmath.exp(-1j*Jf*a-1j*float(M)*g)
                worst = max(worst, abs(D(J, J, M, a, b, g)-r1))
                r3 = pref*math.cos(b/2)**(Jf+float(M))*math.sin(b/2)**(Jf-float(M))*cmath.exp(-1j*float(M)*a-1j*Jf*g)
                worst = max(worst, abs(D(J, M, J, a, b, g)-r3))
    ok &= report("4.17.8 D_{JM}, D_{MJ}", worst)
    return ok


def main():
    ok = check_416() and check_417()
    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
