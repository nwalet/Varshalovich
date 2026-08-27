#!/usr/bin/env python3
r"""Verify Sec 5.6.2 (expansions of products of Y) vs mpmath.spherharm.
Clebsch-Gordan coefficients from sympy.  Tensor product convention (Sec 3.1):
  {Y_l1 (x) Y_l2}_{LM}(Om) = sum_{m1,m2} C^{LM}_{l1 m1 l2 m2} Y_{l1 m1}(Om) Y_{l2 m2}(Om).
"""
import mpmath as mp
from sympy.physics.wigner import clebsch_gordan
mp.mp.dps = 25
pi = mp.pi
_cg = {}
def CG(j1,m1,j2,m2,j3,m3):
    k = (j1,m1,j2,m2,j3,m3)
    if k not in _cg: _cg[k] = mp.mpf(float(clebsch_gordan(j1,j2,j3,m1,m2,m3)))
    return _cg[k]
def Y(l,m,th,ph):
    return mp.spherharm(l,m,th,ph) if abs(m)<=l else mp.mpc(0)
def fac(n): return mp.factorial(int(round(n)))
def dfac(n):
    n=int(round(n))
    if n<=0: return mp.mpf(1)
    r=mp.mpf(1)
    while n>1: r*=n; n-=2
    return r
def report(tag,w,tol=mp.mpf('1e-12')):
    ok=w<tol; print(f"  {tag:32s} {'PASS' if ok else 'FAIL'}  worst={mp.nstr(w,3)}"); return ok

# sample directions
OM = [(mp.mpf('0.7'),mp.mpf('1.1')), (mp.mpf('1.9'),mp.mpf('4.2')), (mp.mpf('2.4'),mp.mpf('0.3'))]
ok = True

# 5.6.9  Y_{l1m1} Y_{l2m2} = sum_L sqrt((2l1+1)(2l2+1)/(4pi(2L+1)))
#         C^{L0}_{l1 0 l2 0} C^{L,m1+m2}_{l1m1 l2m2} Y_{L,m1+m2}
w = mp.mpf(0)
for l1 in range(0,4):
    for l2 in range(0,4):
        for m1 in range(-l1,l1+1):
            for m2 in range(-l2,l2+1):
                M = m1+m2
                for (th,ph) in OM:
                    s = mp.mpc(0)
                    for L in range(abs(l1-l2), l1+l2+1):
                        if abs(M) > L: continue
                        s += (mp.sqrt((2*l1+1)*(2*l2+1)/(4*pi*(2*L+1)))
                              * CG(l1,0,l2,0,L,0)*CG(l1,m1,l2,m2,L,M)*Y(L,M,th,ph))
                    w = max(w, abs(Y(l1,m1,th,ph)*Y(l2,m2,th,ph) - s))
ok &= report("5.6.9 CG series", w)

# 5.6.10  sqrt(4pi(2L+1)/((2l1+1)(2l2+1))) sum_{m1m2} C^{LM}_{l1m1 l2m2} Y_{l1m1}Y_{l2m2}
#          = C^{L0}_{l1 0 l2 0} Y_LM
w = mp.mpf(0)
for l1 in range(0,4):
    for l2 in range(0,4):
        for L in range(abs(l1-l2), l1+l2+1):
            for M in range(-L,L+1):
                for (th,ph) in OM:
                    s = mp.mpc(0)
                    for m1 in range(-l1,l1+1):
                        m2 = M-m1
                        if abs(m2) > l2: continue
                        s += CG(l1,m1,l2,m2,L,M)*Y(l1,m1,th,ph)*Y(l2,m2,th,ph)
                    lhs = mp.sqrt(4*pi*(2*L+1)/((2*l1+1)*(2*l2+1)))*s
                    w = max(w, abs(lhs - CG(l1,0,l2,0,L,0)*Y(L,M,th,ph)))
ok &= report("5.6.10 inverse relation", w)

# 5.6.11  Y1 Y2 Y3 = sum_{L M L' M'} sqrt((2l1+1)(2l2+1)(2l3+1)/((4pi)^2(2L+1)))
#          C^{L'0}_{l1 0 l2 0} C^{L0}_{L' 0 l3 0} C^{L'M'}_{l1m1 l2m2} C^{LM}_{L'M' l3m3} Y_LM
w = mp.mpf(0)
for l1 in range(0,3):
    for l2 in range(0,3):
        for l3 in range(0,3):
            for m1 in range(-l1,l1+1):
                for m2 in range(-l2,l2+1):
                    for m3 in range(-l3,l3+1):
                        M = m1+m2+m3; Mp = m1+m2
                        for (th,ph) in OM:
                            s = mp.mpc(0)
                            for Lp in range(abs(l1-l2), l1+l2+1):
                                if abs(Mp) > Lp: continue
                                for L in range(abs(Lp-l3), Lp+l3+1):
                                    if abs(M) > L: continue
                                    s += (mp.sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/((4*pi)**2*(2*L+1)))
                                          * CG(l1,0,l2,0,Lp,0)*CG(Lp,0,l3,0,L,0)
                                          * CG(l1,m1,l2,m2,Lp,Mp)*CG(Lp,Mp,l3,m3,L,M)*Y(L,M,th,ph))
                            w = max(w, abs(Y(l1,m1,th,ph)*Y(l2,m2,th,ph)*Y(l3,m3,th,ph) - s))
ok &= report("5.6.11 three-Y product", w)

# 5.6.14  {Y_l1 (x) Y_l2}_LM = sqrt((2l1+1)(2l2+1)/(4pi(2L+1))) C^{L0}_{l1 0 l2 0} Y_LM
def tens2(l1,l2,L,M,th,ph):
    s = mp.mpc(0)
    for m1 in range(-l1,l1+1):
        m2 = M-m1
        if abs(m2) > l2: continue
        s += CG(l1,m1,l2,m2,L,M)*Y(l1,m1,th,ph)*Y(l2,m2,th,ph)
    return s
w = mp.mpf(0)
for l1 in range(0,4):
    for l2 in range(0,4):
        for L in range(abs(l1-l2), l1+l2+1):
            for M in range(-L,L+1):
                for (th,ph) in OM:
                    rhs = mp.sqrt((2*l1+1)*(2*l2+1)/(4*pi*(2*L+1)))*CG(l1,0,l2,0,L,0)*Y(L,M,th,ph)
                    w = max(w, abs(tens2(l1,l2,L,M,th,ph) - rhs))
ok &= report("5.6.14 tensor product (2)", w)

# 5.6.15  {{Y_l1 (x) Y_l2}_L' (x) Y_l3}_LM
#   = sqrt((2l1+1)(2l2+1)(2l3+1)/((4pi)^2(2L+1))) C^{L'0}_{l1 0 l2 0} C^{L0}_{L' 0 l3 0} Y_LM
def tens3(l1,l2,Lp,l3,L,M,th,ph):
    s = mp.mpc(0)
    for Mp in range(-Lp,Lp+1):
        m3 = M-Mp
        if abs(m3) > l3: continue
        s += CG(Lp,Mp,l3,m3,L,M)*tens2(l1,l2,Lp,Mp,th,ph)*Y(l3,m3,th,ph)
    return s
w = mp.mpf(0)
for l1 in range(0,3):
    for l2 in range(0,3):
        for l3 in range(0,3):
            for Lp in range(abs(l1-l2), l1+l2+1):
                for L in range(abs(Lp-l3), Lp+l3+1):
                    for M in range(-L,L+1):
                        for (th,ph) in OM[:2]:
                            rhs = (mp.sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/((4*pi)**2*(2*L+1)))
                                   * CG(l1,0,l2,0,Lp,0)*CG(Lp,0,l3,0,L,0)*Y(L,M,th,ph))
                            w = max(w, abs(tens3(l1,l2,Lp,l3,L,M,th,ph) - rhs))
ok &= report("5.6.15 tensor product (3)", w)

# 5.6.17  iterated {Y_1 x Y_1 x ... }_{n m} = sqrt(3^n/(4pi)^{n-1} n!/(2n+1)!!) Y_nm
def iterated(n, M, th, ph):
    # T_1 = Y_1;  T_k = {T_{k-1} (x) Y_1}_k
    Tprev = {mu: Y(1,mu,th,ph) for mu in (-1,0,1)}   # level 1: L=1
    for k in range(2, n+1):
        Tnew = {}
        for Mk in range(-k, k+1):
            s = mp.mpc(0)
            for mu in (-1,0,1):
                Mprev = Mk-mu
                if abs(Mprev) > k-1: continue
                s += CG(k-1,Mprev,1,mu,k,Mk)*Tprev.get(Mprev,mp.mpc(0))*Y(1,mu,th,ph)
            Tnew[Mk] = s
        Tprev = Tnew
    return Tprev[M]
w = mp.mpf(0)
for n in range(1,6):
    for M in range(-n,n+1):
        for (th,ph) in OM:
            rhs = mp.sqrt(mp.mpf(3)**n/(4*pi)**(n-1)*fac(n)/dfac(2*n+1))*Y(n,M,th,ph)
            w = max(w, abs(iterated(n,M,th,ph) - rhs))
ok &= report("5.6.17 iterated Y_1 coupling", w)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
