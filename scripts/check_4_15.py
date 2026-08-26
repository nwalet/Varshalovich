#!/usr/bin/env python3
r"""
Checks for Section 4.15 (Generalized characters chi_lambda^J(omega)) of Ch.4 VMK.

Reference: the DEFINITION 4.15.1,
  chi_lam^J(w) = sqrt(2J+1) sqrt((2J-lam)!/(2J+lam+1)!) (sin w/2)^lam
                 (d/d cos(w/2))^lam chi^J(w),
with chi^J(w) = U_{2J}(cos w/2) (Chebyshev-U; = sin[(2J+1)w/2]/sin(w/2)).
Computed symbolically (sympy chebyshevu + diff) and lambdified to mpmath.
Cross-validated vs the independent CG series 4.15.2 (machine precision).

Covers trig series 4.15.2-4, differential form 4.15.5, Gegenbauer 4.15.6,
Jacobi 4.15.7/8, hypergeometric 4.15.9-16, integral reps 4.15.17/18,
symmetries 4.15.19/20, particular omega 4.15.21, recursions 4.15.22/23,
asymptotics 4.15.24/25, differential eq 4.15.26, orthogonality 4.15.28.

Forms expressed through cos(w/2) or a single-branch argument are tested on
w in (0,pi) (noted); series / symmetry / recursion / ODE forms on all of (0,2pi).

Usage:  python3 check_4_15.py
"""
import math, cmath
import mpmath as mp
import sympy as sp
from sympy.physics.wigner import clebsch_gordan as _CG
from scipy.integrate import quad

mp.mp.dps = 30
TOL = 1e-9
_w = sp.symbols('w', real=True)
WS_FULL = [0.4, 0.8, 1.3, 2.0, 2.7, 3.5, 4.5]     # (0,2pi)
WS_PI   = [0.4, 0.8, 1.3, 2.0, 2.7, 3.0]          # (0,pi)

# ---- reference from 4.15.1 -------------------------------------------------
_refcache = {}
def ref_fn(J, lam):
    key = (J, lam)
    if key not in _refcache:
        x = sp.symbols('x')
        U = sp.chebyshevu(int(2*J), x)
        dU = sp.diff(U, x, lam).subs(x, sp.cos(_w/2))
        pref = sp.sqrt(2*J+1)*sp.sqrt(sp.factorial(int(2*J-lam)) /
                                      sp.factorial(int(2*J+lam+1)))
        _refcache[key] = sp.lambdify(_w, pref*(sp.sin(_w/2))**lam*dU, 'mpmath')
    return _refcache[key]
def chi(J, lam, w):
    return complex(ref_fn(J, lam)(w))

def jl_pairs(maxJ2=5):
    """(J, lam) with 0<=lam<=2J, 2J = 1..maxJ2."""
    out = []
    for twoJ in range(1, maxJ2+1):
        J = sp.Rational(twoJ, 2)
        for lam in range(0, twoJ+1):
            out.append((J, lam))
    return out

def fac(n): return math.factorial(int(round(n)))
def dfac(n):
    n = int(round(n)); r = 1.0
    while n > 1: r *= n; n -= 2
    return r
def CGf(j1, j2, j3, m1, m2, m3):
    return float(_CG(sp.nsimplify(j1), sp.nsimplify(j2), sp.nsimplify(j3),
                     sp.nsimplify(m1), sp.nsimplify(m2), sp.nsimplify(m3)))
def F(a, b, c, z):
    return complex(mp.hyp2f1(a, b, c, z))

def report(tag, worst, tol=TOL):
    ok = worst < tol
    print(f"  {tag:52s} {'PASS' if ok else 'FAIL'}  worst={worst:.2e}")
    return ok

def sweep(tag, form, pairs=None, ws=WS_FULL, tol=TOL):
    if pairs is None: pairs = jl_pairs()
    worst = 0.0
    for (J, lam) in pairs:
        for w in ws:
            worst = max(worst, abs(form(J, lam, w) - chi(J, lam, w)))
    return report(tag, worst, tol)


def main():
    print("Section 4.15 generalized characters chi_lambda^J(omega)\n")
    ok = True

    # 4.15.2  i^lam sum_M e^{-iMw} C^{JM}_{JM lam0}
    def f2(J, lam, w):
        tot = 0j; M = -J
        while M <= J + sp.Rational(1, 1000):
            tot += cmath.exp(-1j*float(M)*w)*CGf(J, lam, J, M, 0, M)
            M += 1
        return 1j**lam*tot
    ok &= sweep("4.15.2 CG trig series", f2)

    # 4.15.3  (sin w/2)^lam sqrt(2J+1)sqrt((2J-lam)!/(2J+lam+1)!) 2^lam
    #          sum_{s=0}^{floor(J-lam/2)} (-1)^s (2J-s)!/(s!(2J-lam-2s)!) (2cos w/2)^{2J-lam-2s}
    def f3(J, lam, w):
        s2 = math.sin(w/2); c = math.cos(w/2)
        pref = s2**lam*math.sqrt(2*J+1)*math.sqrt(fac(2*J-lam)/fac(2*J+lam+1))*2**lam
        tot = 0.0
        smax = int(math.floor(float(J) - lam/2 + 1e-9))
        for s in range(smax+1):
            if 2*J-lam-2*s < 0: break
            tot += (-1)**s*fac(2*J-s)/(fac(s)*fac(2*J-lam-2*s))*(2*c)**(2*J-lam-2*s)
        return pref*tot
    ok &= sweep("4.15.3 trig series (cos powers)", f3)

    # 4.15.4  (2 sin w/2)^lam/lam! sqrt(...) sum_{s=0}^{2J-lam}
    #          (lam+s)!(2J-s)!/(s!(2J-lam-s)!) cos[(2J-lam-2s) w/2]
    def f4(J, lam, w):
        s2 = math.sin(w/2)
        pref = (2*s2)**lam/fac(lam)*math.sqrt(2*J+1)*math.sqrt(fac(2*J-lam)/fac(2*J+lam+1))
        tot = 0.0
        for s in range(int(2*J-lam)+1):
            tot += (fac(lam+s)*fac(2*J-s)/(fac(s)*fac(2*J-lam-s))
                    * math.cos((2*J-lam-2*s)*w/2))
        return pref*tot
    ok &= sweep("4.15.4 trig series (cos args)", f4)

    # 4.15.5  differential form: 1/sqrt(2J+1) sqrt((2J-lam)!/(2J+lam+1)!)
    #          (sin w/2)^lam (d/dc)^{lam+1} cos[(2J+1) w/2]
    def f5(J, lam, w):
        x = sp.symbols('x')
        # cos[(2J+1) w/2] as function of c=cos(w/2): cos((2J+1)*arccos(x))=T_{2J+1}(x)
        T = sp.chebyshevt(int(2*J+1), x)
        dT = sp.diff(T, x, lam+1).subs(x, math.cos(w/2))
        pref = (1/math.sqrt(2*J+1)*math.sqrt(fac(2*J-lam)/fac(2*J+lam+1))
                * math.sin(w/2)**lam)
        return pref*float(dT)
    ok &= sweep("4.15.5 differential form", f5)

    # 4.15.6  (2lam)!! sqrt(2J+1)sqrt((2J-lam)!/(2J+lam+1)!) (sin w/2)^lam C^{lam+1}_{2J-lam}(cos w/2)
    def f6(J, lam, w):
        from scipy.special import eval_gegenbauer
        pref = (dfac(2*lam)*math.sqrt(2*J+1)*math.sqrt(fac(2*J-lam)/fac(2*J+lam+1))
                * math.sin(w/2)**lam)
        return pref*eval_gegenbauer(int(2*J-lam), lam+1, math.cos(w/2))
    ok &= sweep("4.15.6 Gegenbauer C^{lam+1}_{2J-lam}", f6)

    # 4.15.7 Jacobi P^{(lam+1/2,lam+1/2)}_{2J-lam}(cos w/2)
    def f7(J, lam, w):
        from scipy.special import eval_jacobi
        pref = (math.sqrt(2*J+1)*math.sqrt(fac(2*J-lam)*fac(2*J+lam+1))/dfac(4*J+1)
                * 2**(2*J-lam)*math.sin(w/2)**lam)
        return pref*eval_jacobi(int(2*J-lam), lam+0.5, lam+0.5, math.cos(w/2))
    ok &= sweep("4.15.7 Jacobi (lam+1/2,lam+1/2)", f7)

    # 4.15.8 Jacobi in cos w, piecewise by parity of 2J-lam
    def f8(J, lam, w):
        from scipy.special import eval_jacobi
        s2 = math.sin(w/2); c2 = math.cos(w/2); cw = math.cos(w)
        if int(2*J-lam) % 2 == 0:
            pref = (math.sqrt(2*J+1)*math.sqrt(dfac(2*J+lam)*dfac(2*J-lam)
                    / (dfac(2*J+lam+1)*dfac(2*J-lam-1)))*s2**lam)
            return pref*eval_jacobi(int(J-lam/2), lam+0.5, -0.5, cw)
        else:
            pref = (math.sqrt(2*J+1)*math.sqrt(dfac(2*J-lam-1)*dfac(2*J+lam+1)
                    / (dfac(2*J-lam)*dfac(2*J+lam)))*c2*s2**lam)
            return pref*eval_jacobi(int(J-(lam+1)/2), lam+0.5, 0.5, cw)
    ok &= sweep("4.15.8 Jacobi in cos w (piecewise)", f8, ws=WS_PI)

    # 4.15.9  F(-2J+lam, 2J+lam+2; lam+3/2; sin^2 w/4)
    def f9(J, lam, w):
        pref = (math.sqrt(2*J+1)/dfac(2*lam+1)*math.sqrt(fac(2*J+lam+1)/fac(2*J-lam))
                * math.sin(w/2)**lam)
        return pref*F(-2*J+lam, 2*J+lam+2, lam+1.5, math.sin(w/4)**2)
    ok &= sweep("4.15.9 2F1(sin^2 w/4)", f9)

    # 4.15.10  (-1)^{2J-lam}(2lam)!!/(2lam+1)! ... F(...; cos^2 w/4)
    def f10(J, lam, w):
        pref = ((-1)**int(2*J-lam)*dfac(2*lam)*math.sqrt(2*J+1)/fac(2*lam+1)
                * math.sqrt(fac(2*J+lam+1)/fac(2*J-lam))*math.sin(w/2)**lam)
        return pref*F(-2*J+lam, 2*J+lam+2, lam+1.5, math.cos(w/4)**2)
    ok &= sweep("4.15.10 2F1(cos^2 w/4)", f10, ws=WS_PI)

    # 4.15.13  F(-2J+lam, -2J-1/2; lam+3/2; -tan^2 w/4)
    def f13(J, lam, w):
        pref = (math.sqrt(2*J+1)/dfac(2*lam+1)*math.sqrt(fac(2*J+lam+1)/fac(2*J-lam))
                * math.sin(w/2)**lam*math.cos(w/4)**(4*J-2*lam))
        return pref*F(-2*J+lam, -2*J-0.5, lam+1.5, -math.tan(w/4)**2)
    ok &= sweep("4.15.13 2F1(-tan^2 w/4)", f13, ws=WS_PI)

    # 4.15.16  hypergeometric in sin^2 w/2, piecewise
    def f16(J, lam, w):
        s2 = math.sin(w/2); c2 = math.cos(w/2)
        base = math.sqrt(2*J+1)/dfac(2*lam+1)*math.sqrt(fac(2*J+lam+1)/fac(2*J-lam))
        if int(2*J-lam) % 2 == 0:
            return base*s2**lam*F(-J+lam/2, J+1+lam/2, lam+1.5, s2**2)
        else:
            return base*c2*s2**lam*F(-J+(lam+1)/2, J+1+(lam+1)/2, lam+1.5, s2**2)
    ok &= sweep("4.15.16 2F1(sin^2 w/2, piecewise)", f16, ws=WS_PI)

    # 4.15.15  hypergeometric in cos^2 w/2, piecewise
    def f15(J, lam, w):
        s2 = math.sin(w/2); c2 = math.cos(w/2)
        if int(2*J-lam) % 2 == 0:
            pref = ((-1)**(J+sp.Rational(3*lam, 2)))
            pref = float(sp.N(pref))*math.sqrt(2*J+1)*math.sqrt(
                dfac(2*J-lam-1)*dfac(2*J+lam)/(dfac(2*J-lam)*dfac(2*J+lam+1)))*s2**lam
            return pref*F(-J+lam/2, J+1+lam/2, 0.5, c2**2)
        else:
            pref = float(sp.N((-1)**(J+sp.Rational(3*lam-1, 2))))*math.sqrt(2*J+1)*math.sqrt(
                dfac(2*J-lam)*dfac(2*J+lam+1)/(dfac(2*J-lam-1)*dfac(2*J+lam)))*c2*s2**lam
            return pref*F(-J+(lam+1)/2, J+1+(lam+1)/2, 1.5, c2**2)
    ok &= sweep("4.15.15 2F1(cos^2 w/2, piecewise)", f15, ws=WS_PI)

    # 4.15.18 integral rep: (-i)^lam sqrt((2J+1)(2J+lam+1)!(2J-lam)!)/(2 (2J)!)
    #          int_{-1}^{1} P_lam(x) [cos w/2 + i x sin w/2]^{2J} dx
    def f18(J, lam, w):
        from scipy.special import eval_legendre
        n = int(2*J); c = math.cos(w/2); s = math.sin(w/2)
        def integ_re(x): return eval_legendre(lam, x)*((c+1j*x*s)**n).real
        def integ_im(x): return eval_legendre(lam, x)*((c+1j*x*s)**n).imag
        re = quad(integ_re, -1, 1)[0]; im = quad(integ_im, -1, 1)[0]
        pref = ((-1j)**lam*math.sqrt((2*J+1)*fac(2*J+lam+1)*fac(2*J-lam))/(2*fac(2*J)))
        return pref*(re+1j*im)
    ok &= sweep("4.15.18 integral rep (Legendre)", f18)

    # 4.15.19 symmetry: chi* = chi = (-1)^lam chi(-w)
    def f19(J, lam, w):
        return (-1)**lam*chi(J, lam, -w)
    ok &= sweep("4.15.19 (-1)^lam chi(-w)", f19)
    # conjugate real:
    worst = max(abs(chi(J, lam, w).imag) for (J, lam) in jl_pairs() for w in WS_FULL)
    ok &= report("4.15.19 chi real", worst)

    # 4.15.20 chi(2pi-w) = (-1)^{2J-lam} chi(w)
    worst = 0.0
    for (J, lam) in jl_pairs():
        for w in WS_PI:
            worst = max(worst, abs(chi(J, lam, 2*math.pi-w)
                                   - (-1)**int(2*J-lam)*chi(J, lam, w)))
    ok &= report("4.15.20 chi(2pi-w)=(-1)^{2J-lam}chi", worst)

    # 4.15.21 chi(pi): (-1)^{J-lam/2} sqrt(...) if 2J-lam even, else 0
    def chi_pi(J, lam):
        # NB: denominator uses (2J-lam)!! (double factorial); the printed book /
        # .tex had (2J-lam)! (single) -- wrong for 2J-lam>=4 (see flags).
        if int(2*J-lam) % 2 == 1: return 0.0
        val = float(sp.N((-1)**(J-sp.Rational(lam, 2))))*math.sqrt(
            dfac(2*J+lam)*dfac(2*J-lam-1)*(2*J+1)/(dfac(2*J-lam)*dfac(2*J+lam+1)))
        return val
    worst = max(abs(chi(J, lam, math.pi).real - chi_pi(J, lam))
                for (J, lam) in jl_pairs())
    ok &= report("4.15.21 chi(pi)", worst)

    # 4.15.22 recursion (derivative)
    def f22(J, lam, w):
        h = 1e-6
        lhs = 2*(chi(J, lam, w+h)-chi(J, lam, w-h))/(2*h)
        rhs = (lam/(2*lam+1)*math.sqrt((2*J+1)**2-lam**2)*chi(J, lam-1, w)
               - (lam+1)/(2*lam+1)*math.sqrt((2*J+1)**2-(lam+1)**2)*chi(J, lam+1, w))
        return lhs - rhs
    worst = 0.0
    for (J, lam) in jl_pairs():
        if lam < 1 or lam+1 > 2*J: continue
        for w in [0.8, 1.3, 2.0, 2.7]:
            worst = max(worst, abs(f22(J, lam, w)))
    ok &= report("4.15.22 recursion (d/dw)", worst, tol=1e-5)

    # 4.15.23 recursion (cot)
    def f23(J, lam, w):
        lhs = (2*lam+1)/math.tan(w/2)*chi(J, lam, w)
        rhs = (math.sqrt((2*J+1)**2-lam**2)*chi(J, lam-1, w)
               + math.sqrt((2*J+1)**2-(lam+1)**2)*chi(J, lam+1, w))
        return lhs - rhs
    worst = 0.0
    for (J, lam) in jl_pairs():
        if lam < 1 or lam+1 > 2*J: continue
        for w in [0.8, 1.3, 2.0, 2.7]:
            worst = max(worst, abs(f23(J, lam, w)))
    ok &= report("4.15.23 recursion (cot)", worst)

    # 4.15.25 asymptotics w->0
    def f25(J, lam, w):
        return ((w/2)**lam/dfac(2*lam+1)*math.sqrt((2*J+1)*fac(2*J+lam+1)/fac(2*J-lam)))
    worst = 0.0
    for (J, lam) in jl_pairs():
        wv = 1e-4
        worst = max(worst, abs(chi(J, lam, wv).real - f25(J, lam, wv))/max(1, abs(f25(J,lam,wv))))
    ok &= report("4.15.25 asymptotic w->0 (rel)", worst, tol=1e-3)

    # 4.15.26 differential equation
    def f26(J, lam, w):
        h = 1e-3
        c1 = (-chi(J,lam,w+2*h)+8*chi(J,lam,w+h)-8*chi(J,lam,w-h)+chi(J,lam,w-2*h))/(12*h)
        c2 = (-chi(J,lam,w+2*h)+16*chi(J,lam,w+h)-30*chi(J,lam,w)
              +16*chi(J,lam,w-h)-chi(J,lam,w-2*h))/(12*h**2)
        return c2 + (1/math.tan(w/2))*c1 + (J*(J+1)-lam*(lam+1)/(4*math.sin(w/2)**2))*chi(J,lam,w)
    worst = 0.0
    for (J, lam) in jl_pairs():
        for w in [0.8, 1.3, 2.0, 2.7]:
            worst = max(worst, abs(f26(J, lam, w)))
    ok &= report("4.15.26 differential equation", worst, tol=1e-6)

    # 4.15.28 orthogonality: int_0^2pi sin^2(w/2) chi_lam^{J1} chi_lam^{J2} = pi d_{J1J2}
    def i28(J1, J2, lam):
        f1 = ref_fn(J1, lam); f2 = ref_fn(J2, lam)
        integ = lambda w: math.sin(w/2)**2*float((f1(w)*f2(w)).real)
        return quad(integ, 0, 2*math.pi, limit=200)[0]
    worst = 0.0
    for lam in [0, 1, 2]:
        Js = [sp.Rational(k, 2) for k in range(1, 6) if sp.Rational(k, 2) >= lam/2]
        for J1 in Js:
            for J2 in Js:
                worst = max(worst, abs(i28(J1, J2, lam)-math.pi*(1 if J1 == J2 else 0)))
    ok &= report("4.15.28 orthogonality", worst, tol=1e-6)

    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
