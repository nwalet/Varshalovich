#!/usr/bin/env python3
r"""Verify Sec 5.17 bipolar-harmonic expansions (5.17.20, 5.17.34-5.17.39) vs
mpmath.spherharm.  These expand r^N Y_LM(Omega) (Omega = direction of
r = r1 - r2) in bipolar harmonics {Y_l1(O1) x Y_l2(O2)}_LM.  CG from sympy.
"""
import math, cmath
import numpy as np
import mpmath as mp
from sympy.physics.wigner import clebsch_gordan
mp.mp.dps = 20
pi = mp.pi
def Y(l,m,O): return complex(mp.spherharm(l,m,O[0],O[1])) if abs(m)<=l else 0j
_cg={}
def CG(j1,m1,j2,m2,j3,m3):
    k=(j1,m1,j2,m2,j3,m3)
    if k not in _cg: _cg[k]=float(clebsch_gordan(j1,j2,j3,m1,m2,m3))
    return _cg[k]
def fac(n): return mp.factorial(int(round(n)))
def jl(l,x): return complex(mp.sqrt(pi/(2*x))*mp.besselj(l+mp.mpf(1)/2,x))
def nl(l,x): return complex(mp.sqrt(pi/(2*x))*mp.bessely(l+mp.mpf(1)/2,x))
def Bip(l1,l2,L,M,O1,O2):
    s=0j
    for m1 in range(-l1,l1+1):
        m2=M-m1
        if abs(m2)>l2: continue
        s+=CG(l1,m1,l2,m2,L,M)*Y(l1,m1,O1)*Y(l2,m2,O2)
    return s
def report(tag,w,tol=1e-10):
    w=float(w); ok=w<tol; print(f"  {tag:36s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok
ok=True

def dvec(r,th,ph): return np.array([r*math.sin(th)*math.cos(ph), r*math.sin(th)*math.sin(ph), r*math.cos(th)])
def vdir(v):
    r=np.linalg.norm(v); return r, math.acos(max(-1,min(1,v[2]/r))), math.atan2(v[1],v[0])%(2*math.pi)

# geometry: r1<r2 for the multipole-type expansions
GEOS=[(0.6,0.7,1.1, 1.9,1.9,4.2), (0.5,1.2,0.3, 1.7,0.6,2.0)]  # (r1,th1,ph1, r2,th2,ph2)
LMAX=60

print("Sec 5.17.6/7 bipolar expansions")
# 5.17.35  r^L Y_LM = sqrt(4pi(2L+1)!) sum_{l1+l2=L} (-1)^l2 r1^l1 r2^l2/sqrt((2l1+1)!(2l2+1)!){Y_l1 x Y_l2}_LM
w=0.0
for (r1,t1,p1,r2,t2,p2) in GEOS:
    O1=(t1,p1); O2=(t2,p2)
    rvec=dvec(r1,t1,p1)-dvec(r2,t2,p2); r,th,ph=vdir(rvec)
    for L in range(0,4):
        for M in range(-L,L+1):
            lhs=r**L*Y(L,M,(th,ph))
            s=mp.sqrt(4*pi*fac(2*L+1))*sum((-1)**l2*mp.mpf(r1)**(L-l2)*mp.mpf(r2)**l2
                     /mp.sqrt(fac(2*(L-l2)+1)*fac(2*l2+1))*Bip(L-l2,l2,L,M,O1,O2) for l2 in range(0,L+1))
            w=max(w,abs(lhs-complex(s)))
ok&=report("5.17.35 r^L Y_LM (finite)", w)

# 5.17.36  (1/r^{L+1})Y_LM = sqrt(4pi/(2L)!) sum_{l2-l1=L} (-1)^l2 sqrt((2l2)!/(2l1+1)!) r1^l1/r2^{l2+1}{Y_l1 x Y_l2}_LM
w=0.0
for (r1,t1,p1,r2,t2,p2) in GEOS:
    O1=(t1,p1); O2=(t2,p2)
    rvec=dvec(r1,t1,p1)-dvec(r2,t2,p2); r,th,ph=vdir(rvec)
    for L in range(0,4):
        for M in range(-L,L+1):
            lhs=Y(L,M,(th,ph))/r**(L+1)
            s=mp.sqrt(4*pi/fac(2*L))*sum((-1)**(l1+L)*mp.sqrt(fac(2*(l1+L))/fac(2*l1+1))
                     *mp.mpf(r1)**l1/mp.mpf(r2)**(l1+L+1)*Bip(l1,l1+L,L,M,O1,O2) for l1 in range(0,LMAX))
            w=max(w,abs(lhs-complex(s)))
ok&=report("5.17.36 Y_LM/r^{L+1} (r1<r2)", w)

# 5.17.37  r_mu/r^3 = 4pi sum_{l>=1} (-1)^l sqrt(l/3) r1^{l-1}/r2^{l+1}{Y_{l-1} x Y_l}_{1 mu}, r1<r2
# spherical component r_mu = sqrt(4pi/3) r Y_{1 mu}(th,ph)
w=0.0
for (r1,t1,p1,r2,t2,p2) in GEOS:
    O1=(t1,p1); O2=(t2,p2)
    rvec=dvec(r1,t1,p1)-dvec(r2,t2,p2); r,th,ph=vdir(rvec)
    for mu in (-1,0,1):
        rmu=math.sqrt(4*math.pi/3)*r*Y(1,mu,(th,ph))
        lhs=rmu/r**3
        s=4*pi*sum((-1)**l*mp.sqrt(mp.mpf(l)/3)*mp.mpf(r1)**(l-1)/mp.mpf(r2)**(l+1)
                   *Bip(l-1,l,1,mu,O1,O2) for l in range(1,LMAX))
        w=max(w,abs(lhs-complex(s)))
ok&=report("5.17.37 r_mu/r^3 (r1<r2)", w)

# 5.17.34  z_L(kr)Y_LM = sqrt(4pi/(2L+1)) sum_{l1,l2} i^{l1-l2-L} sqrt((2l1+1)(2l2+1))
#          C^{L0}_{l1 0 l2 0} j_{l1}(k r1) z_{l2}(k r2){Y_l1 x Y_l2}_LM, r1<r2  (test z=n Neumann)
w=0.0
k=1.2
for (r1,t1,p1,r2,t2,p2) in GEOS:
    O1=(t1,p1); O2=(t2,p2)
    rvec=dvec(r1,t1,p1)-dvec(r2,t2,p2); r,th,ph=vdir(rvec)
    for L in range(0,3):
        for M in range(-L,L+1):
            lhs=nl(L,k*r)*Y(L,M,(th,ph))
            s=0j
            for l1 in range(0,40):
                for l2 in range(abs(l1-L),l1+L+1):
                    if (l1+l2+L)%2: continue
                    s+=(1j)**(l1-l2-L)*math.sqrt((2*l1+1)*(2*l2+1))*CG(l1,0,l2,0,L,0)*jl(l1,k*r1)*nl(l2,k*r2)*Bip(l1,l2,L,M,O1,O2)
            s*=math.sqrt(4*math.pi/(2*L+1))
            w=max(w,abs(lhs-s))
ok&=report("5.17.34 spherical wave (z=n, r1<r2)", w)

# 5.17.38/39  r^N Y_LM = 4pi sum_{l1 l2} a^{NL}_{l1l2} sqrt((2l1+1)(2l2+1)/(4pi(2L+1))) C^{L0}_{l1 0 l2 0}{Y_l1 x Y_l2}_LM
# NOTE: 5.17.39's coefficient is a GENERIC-N closed form; it degenerates at
# integer-N pole configurations of its Gamma factors:
#   (L-N)/2 a non-positive integer  [Gamma((L-N)/2) denominator pole] -> the
#     polynomial cases r^N Y_LM (N>=L, N==L mod 2), covered by 5.17.35;
#   (L+N+3)/2 a non-positive integer [Gamma((L+N+3)/2) numerator pole] -> very
#     negative N with L+N+3<=0 (e.g. N=-3,L=0).
# Where both Gammas are regular the formula is exact.  We test such (N,L) pairs.
def aNL(l1,l2,N,L,r1,r2):
    from mpmath import gamma, rf
    pref=((-1)**((l1-l2-L)//2)*mp.mpf(2)**l1/dfacmp(2*l1+1)
          *gamma(mp.mpf(l1+l2-N)/2)*gamma(mp.mpf(L+N+3)/2)
          /(gamma(mp.mpf(L-N)/2)*gamma(mp.mpf(-l1+l2+N+3)/2)))
    return pref*mp.mpf(r2)**N*(mp.mpf(r1)/r2)**l1*mp.hyp2f1(mp.mpf(l1+l2-N)/2, mp.mpf(l1-l2-N-1)/2, l1+mp.mpf(3)/2, mp.mpf(r1)**2/r2**2)
def dfacmp(n):
    n=int(round(n))
    if n<=0: return mp.mpf(1)
    r=mp.mpf(1)
    while n>1: r*=n; n-=2
    return r
def regular(N,L):  # both Gamma args away from non-positive-integer poles
    a=(L-N)/2; b=(L+N+3)/2
    bad=lambda z: (float(z).is_integer() and z<=0)
    return not (bad(a) or bad(b))
NL_PAIRS=[(N,L) for N in [-1,-3,-5,1,3] for L in range(0,3) if regular(N,L)]
w=0.0
for (r1,t1,p1,r2,t2,p2) in GEOS:
    O1=(t1,p1); O2=(t2,p2)
    rvec=dvec(r1,t1,p1)-dvec(r2,t2,p2); r,th,ph=vdir(rvec)
    for (N,L) in NL_PAIRS:
        for M in range(-L,L+1):
            lhs=r**N*Y(L,M,(th,ph))
            s=mp.mpc(0)
            for l1 in range(0,40):
                for l2 in range(abs(l1-L),l1+L+1):
                    if (l1+l2-L)%2: continue
                    try: a=aNL(l1,l2,N,L,r1,r2)
                    except Exception: continue
                    if not mp.isfinite(a.real): continue
                    s+=a*mp.sqrt((2*l1+1)*(2*l2+1)/(4*pi*(2*L+1)))*CG(l1,0,l2,0,L,0)*Bip(l1,l2,L,M,O1,O2)
            w=max(w,abs(lhs-complex(4*pi*s)))
ok&=report(f"5.17.38/39 r^N Y_LM ({len(NL_PAIRS)} regular N,L)", w)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
