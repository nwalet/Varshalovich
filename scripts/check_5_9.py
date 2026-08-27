#!/usr/bin/env python3
r"""Verify Sec 5.9 (integrals involving spherical harmonics) vs
mpmath.spherharm.  Each Y_{lm}(th,ph)=A_{lm}(th) e^{i m ph} with
A_{lm}(th)=spherharm(l,m,th,0); the ph-integral is done analytically
(=> 2pi delta on sum of m's), leaving fast 1-D th-quadratures.
Sec 5.9.2 (5.9.6-5.9.9) are distributional/operator identities (delta(q-k),
L-hat, curl) -- not verified numerically here; notation-fixed only.
"""
import math
import mpmath as mp
from sympy.physics.wigner import wigner_3j, clebsch_gordan
mp.mp.dps = 25
pi = mp.pi
def A(l, m, th):
    if abs(m) > l: return mp.mpf(0)
    return mp.spherharm(l, m, th, 0)          # real: N_lm P_l^m(cos th)
def fac(n): return mp.factorial(int(round(n)))
def dfac(n):
    n = int(round(n))
    if n in (0, -1): return mp.mpf(1)
    if n > 0:
        r = mp.mpf(1)
        while n > 1: r *= n; n -= 2
        return r
    # n negative
    if n % 2 == 0: return mp.inf                 # negative even: pole
    k = (-n - 1)//2                              # n = -(2k+1)
    return (-1)**k / dfac(2*k-1)
def C3j(l1,l2,l3,m1,m2,m3): return mp.mpf(float(wigner_3j(l1,l2,l3,m1,m2,m3)))
def CG(l1,m1,l2,m2,l3,m3): return mp.mpf(float(clebsch_gordan(l1,l2,l3,m1,m2,m3)))
def th_int(f): return mp.quad(f, [0, pi])
def report(tag, w, tol=mp.mpf('1e-12')):
    ok = w < tol; print(f"  {tag:36s} {'PASS' if ok else 'FAIL'}  worst={mp.nstr(w,3)}"); return ok

ok = True
print("Sec 5.9.1 solid-angle integrals")
# 5.9.1  int int sin Y_lm = sqrt(4pi) d_{l0} d_{m0}
w = mp.mpf(0)
for l in range(0,5):
    for m in range(-l,l+1):
        val = (2*pi if m==0 else 0)*th_int(lambda th,l=l: mp.sin(th)*A(l,0,th)) if m==0 else mp.mpf(0)
        rhs = mp.sqrt(4*pi) if (l==0 and m==0) else mp.mpf(0)
        w = max(w, abs(val-rhs))
ok &= report("5.9.1 int Y", w)
# 5.9.2  int int sin Y_{l1m1} Y*_{l2m2} = d_{l1l2} d_{m1m2}
w = mp.mpf(0)
for l1 in range(0,4):
    for l2 in range(0,4):
        for m in range(-min(l1,l2), min(l1,l2)+1):
            val = 2*pi*th_int(lambda th,l1=l1,l2=l2,m=m: mp.sin(th)*A(l1,m,th)*A(l2,m,th))
            rhs = mp.mpf(1) if l1==l2 else mp.mpf(0)
            w = max(w, abs(val-rhs))
ok &= report("5.9.2 orthonormality", w)
# 5.9.3  int int sin Y_{l1m1} Y_{l2m2} = (-1)^{m2} d_{l1l2} d_{-m1 m2}
w = mp.mpf(0)
for l1 in range(0,4):
    for l2 in range(0,4):
        for m1 in range(-l1,l1+1):
            m2 = -m1
            if abs(m2) > l2: continue
            val = 2*pi*th_int(lambda th,l1=l1,l2=l2,m1=m1,m2=m2: mp.sin(th)*A(l1,m1,th)*A(l2,m2,th))
            rhs = (-1)**m2 if l1==l2 else mp.mpf(0)
            w = max(w, abs(val-rhs))
ok &= report("5.9.3 int Y Y", w)
# 5.9.4  int int sin Y_{l1m1}Y_{l2m2}Y*_{l3m3}
#  = sqrt((2l1+1)(2l2+1)/(4pi(2l3+1))) C^{l3 0}_{l1 0 l2 0} C^{l3 m3}_{l1 m1 l2 m2}
w = mp.mpf(0)
for l1 in range(0,3):
    for l2 in range(0,3):
        for l3 in range(abs(l1-l2), l1+l2+1):
            for m1 in range(-l1,l1+1):
                for m2 in range(-l2,l2+1):
                    m3 = m1+m2
                    if abs(m3) > l3: continue
                    val = 2*pi*th_int(lambda th,a=l1,b=l2,c=l3,p=m1,q=m2,r=m3:
                                      mp.sin(th)*A(a,p,th)*A(b,q,th)*A(c,r,th))
                    rhs = (mp.sqrt((2*l1+1)*(2*l2+1)/(4*pi*(2*l3+1)))
                           * CG(l1,0,l2,0,l3,0)*CG(l1,m1,l2,m2,l3,m3))
                    w = max(w, abs(val-rhs))
ok &= report("5.9.4 three-Y (CG)", w)
# 5.9.5  int int sin Y_{l1m1}Y_{l2m2}Y_{l3m3}
#  = sqrt((2l1+1)(2l2+1)(2l3+1)/4pi) (l1 l2 l3;0 0 0)(l1 l2 l3;m1 m2 m3)
w = mp.mpf(0)
for l1 in range(0,3):
    for l2 in range(0,3):
        for l3 in range(abs(l1-l2), l1+l2+1):
            for m1 in range(-l1,l1+1):
                for m2 in range(-l2,l2+1):
                    m3 = -(m1+m2)
                    if abs(m3) > l3: continue
                    val = 2*pi*th_int(lambda th,a=l1,b=l2,c=l3,p=m1,q=m2,r=m3:
                                      mp.sin(th)*A(a,p,th)*A(b,q,th)*A(c,r,th))
                    rhs = (mp.sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*pi))
                           * C3j(l1,l2,l3,0,0,0)*C3j(l1,l2,l3,m1,m2,m3))
                    w = max(w, abs(val-rhs))
ok &= report("5.9.5 three-Y (3j)", w)

print("\nSec 5.9.3 theta-integrals")
# 5.9.10  int_0^pi A_lm A_l'm sin dth = d_{ll'}/(2pi)
w = mp.mpf(0)
for m in range(0,3):
    for l in range(m,5):
        for lp in range(m,5):
            val = th_int(lambda th,l=l,lp=lp,m=m: A(l,m,th)*A(lp,m,th)*mp.sin(th))
            rhs = (1/(2*pi)) if l==lp else mp.mpf(0)
            w = max(w, abs(val-rhs))
ok &= report("5.9.10 orthogonality (sin)", w)
# 5.9.11  int_0^pi A_lm A_lm' dth/sin = (2l+1)/(4pi m) d_{mm'} (m,m'>0)
w = mp.mpf(0)
for l in range(1,5):
    for m in range(1,l+1):
        for mp_ in range(1,l+1):
            val = th_int(lambda th,l=l,m=m,mp_=mp_: A(l,m,th)*A(l,mp_,th)/mp.sin(th))
            rhs = ((2*l+1)/(4*pi*m)) if m==mp_ else mp.mpf(0)
            w = max(w, abs(val-rhs))
ok &= report("5.9.11 orthogonality (1/sin)", w)
# 5.9.12  int_0^{pi/2} sin^{m+1} cos^n A_lm dth
#  = sqrt((2l+1)/4pi (l+m)!/(l-m)!) (-1)^m n!/((n+l+m+1)!!(n-l+m)!!)
w = mp.mpf(0)
for l in range(0,6):
    for m in range(0,l+1):
        for n in range(0,6):
            val = mp.quad(lambda th,l=l,m=m,n=n: mp.sin(th)**(m+1)*mp.cos(th)**n*A(l,m,th),
                          [0, pi/2])
            dd = dfac(n+l+m+1)*dfac(n-l+m)
            rhs = (mp.mpf(0) if (dd==mp.inf or dd==-mp.inf)
                   else mp.sqrt((2*l+1)/(4*pi)*fac(l+m)/fac(l-m))*(-1)**m*fac(n)/dd)
            w = max(w, abs(val-rhs))
ok &= report("5.9.12 sin^{m+1}cos^n integral", w, tol=mp.mpf('1e-10'))
# 5.9.13  Wronskian identity: LHS int over d(cos th) = RHS boundary term.
# BOOK MISPRINT: printed closed-form coefficients (l2+m2),(l1+m1) are the
# un-normalized P_l^m form; for VMK's normalized Y_lm they must be
# K(l,m)=sqrt((2l+1)(l^2-m^2)/(2l-1)).  Verified for general m1,m2.
def Kc(l, m): return mp.sqrt((2*l+1)*(l*l-m*m)/(2*l-1))
def integrand(x, l1,m1,l2,m2):
    th = mp.acos(x)
    return ((l1-l2)*(l1+l2+1) - (m1*m1-m2*m2)/(1-x*x))*A(l1,m1,th)*A(l2,m2,th)
def bdry(x, l1,m1,l2,m2):
    th = mp.acos(x)
    return (x*(l1-l2)*A(l1,m1,th)*A(l2,m2,th)
            + Kc(l2,m2)*A(l1,m1,th)*A(l2-1,m2,th)
            - Kc(l1,m1)*A(l1-1,m1,th)*A(l2,m2,th))
w = mp.mpf(0)
xa, xb = mp.mpf('-0.6'), mp.mpf('0.7')
for l1 in range(1,4):
    for l2 in range(1,4):
        for m1 in range(0, l1+1):
            for m2 in range(0, l2+1):
                lhs = mp.quad(lambda x,l1=l1,m1=m1,l2=l2,m2=m2: integrand(x,l1,m1,l2,m2), [xa, xb])
                rhs = bdry(xb,l1,m1,l2,m2) - bdry(xa,l1,m1,l2,m2)
                w = max(w, abs(lhs-rhs))
ok &= report("5.9.13 Wronskian [corr K coeff]", w, tol=mp.mpf('1e-10'))

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
