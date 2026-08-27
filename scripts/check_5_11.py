#!/usr/bin/env python3
r"""Verify Sec 5.11 (generating functions for Y_lm) vs mpmath.spherharm.
Term-wise convergent generating functions 5.11.1/2/3 are summed directly
(|t|<1 or |t|>1).  The step-function relations 5.11.4/5/7 have LHS
singularities at vartheta=psi and are only distributionally (Abel) summable;
5.11.7 (m=0, Legendre) is conditionally convergent and confirmed here by Abel
summation.  5.11.4/5 are the m>0 analogues (scan-matched, PDF p.165).
"""
import mpmath as mp
import numpy as np
mp.mp.dps = 30
pi = mp.pi
def Y(l, m, th, ph):
    return mp.spherharm(l, m, th, ph) if abs(m) <= l else mp.mpc(0)
def fac(n): return mp.factorial(int(round(n)))
def dfac(n):
    n = int(round(n))
    if n <= 0: return mp.mpf(1)
    r = mp.mpf(1)
    while n > 1: r *= n; n -= 2
    return r
def P(l, x): return mp.legendre(l, x)
def report(tag, w, tol=mp.mpf('1e-10')):
    ok = w < tol; print(f"  {tag:34s} {'PASS' if ok else 'FAIL'}  worst={mp.nstr(w,3)}"); return ok

LMAX = 130
TH = [mp.mpf('0.4'), mp.mpf('1.3'), mp.mpf('2.2')]
ok = True

# 5.11.1  1/R^{2m+1} = (-1)^m/((2m-1)!!(sin)^m) sum ... Y_lm(th,0)
#   |t|<1: t^{l-m};  |t|>1: 1/t^{l+m+1}
def s111(m, t, th, big):
    s = mp.mpc(0)
    for l in range(m, LMAX+1):
        wgt = (1/t**(l+m+1)) if big else t**(l-m)
        s += wgt*mp.sqrt(4*pi/(2*l+1)*fac(l+m)/fac(l-m))*Y(l,m,th,0)
    return (-1)**m/(dfac(2*m-1)*mp.sin(th)**m)*s
w = mp.mpf(0)
for m in range(0,4):
    for t in (mp.mpf('0.3'), mp.mpf('0.6')):
        for th in TH:
            R = mp.sqrt(1-2*t*mp.cos(th)+t*t)
            w = max(w, abs(s111(m,t,th,False) - 1/R**(2*m+1)))
    for t in (mp.mpf('1.8'), mp.mpf('3.0')):
        for th in TH:
            R = mp.sqrt(1-2*t*mp.cos(th)+t*t)
            w = max(w, abs(s111(m,t,th,True) - 1/R**(2*m+1)))
ok &= report("5.11.1 1/R^{2m+1} (both |t|)", w)

# 5.11.2  1/R = sum t^l P_l (|t|<1);  sum t^{-(l+1)} P_l (|t|>1)
w = mp.mpf(0)
for t in (mp.mpf('0.3'), mp.mpf('0.7')):
    for th in TH:
        s = sum(t**l*P(l,mp.cos(th)) for l in range(0,LMAX+1))
        w = max(w, abs(s - 1/mp.sqrt(1-2*t*mp.cos(th)+t*t)))
for t in (mp.mpf('1.8'), mp.mpf('3.0')):
    for th in TH:
        s = sum(1/t**(l+1)*P(l,mp.cos(th)) for l in range(0,LMAX+1))
        w = max(w, abs(s - 1/mp.sqrt(1-2*t*mp.cos(th)+t*t)))
ok &= report("5.11.2 1/R, m=0 (both |t|)", w)

# 5.11.3  [(1+R)^2-t^2]^{-m}/R = (-1)^m/(2^m (sin)^m) sum t^{l-m} l! sqrt(4pi/((2l+1)(l+m)!(l-m)!)) Y
w = mp.mpf(0)
for m in range(0,4):
    for t in (mp.mpf('0.3'), mp.mpf('0.6')):
        for th in TH:
            s = mp.mpc(0)
            for l in range(m, LMAX+1):
                s += t**(l-m)*fac(l)*mp.sqrt(4*pi/((2*l+1)*fac(l+m)*fac(l-m)))*Y(l,m,th,0)
            R = mp.sqrt(1-2*t*mp.cos(th)+t*t)
            lhs = ((1+R)**2 - t*t)**(-m)/R
            w = max(w, abs(lhs - (-1)**m/(2**m*mp.sin(th)**m)*s))
ok &= report("5.11.3 [(1+R)^2-t^2]^{-m}/R", w)

# 5.11.7  distributional; Abel-summed RHS -> LHS as r->1 (numpy, m=0 Legendre)
def Pnp(x, L):
    Pl = np.zeros(L+1); Pl[0] = 1.0
    if L >= 1: Pl[1] = x
    for l in range(1, L):
        Pl[l+1] = ((2*l+1)*x*Pl[l] - l*Pl[l-1])/(l+1)
    return Pl
worst = 0.0
for (th, ps) in ((1.0,0.6), (0.7,1.5), (2.0,1.2)):
    L = 20000; x = np.cos(th); Pl = Pnp(x, L); l = np.arange(L+1)
    val = np.sqrt(2)*(0.9997**l*np.exp(1j*(2*l+1)*ps/2)*Pl).sum()
    cph, cth = np.cos(ps), np.cos(th)
    a = 1/np.sqrt(cph-cth) if cph >= cth else 0.0
    b = 1/np.sqrt(cth-cph) if cth >= cph else 0.0
    lhs = a + 1j*b
    worst = max(worst, abs(val-lhs))
print(f"  {'5.11.7 Abel-sum (r=0.9997)':34s} {'PASS' if worst<3e-3 else 'FAIL'}  worst={worst:.2e}")
ok &= worst < 3e-3
print("  5.11.4/5  distributional (Theta + singular LHS): scan-matched, not term-wise convergent")

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
