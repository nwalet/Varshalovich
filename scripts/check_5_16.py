#!/usr/bin/env python3
r"""Verify Sec 5.16 (bipolar and tripolar spherical harmonics) vs
mpmath.spherharm.  CG/3j/6j/9j from sympy.  Bipolar:
  {Y_l1(O1) x Y_l2(O2)}_{LM} = sum_{m1 m2} C^{LM}_{l1 m1 l2 m2} Y_{l1 m1}(O1) Y_{l2 m2}(O2).
Tripolar (l1 x (l2 x l3)_lam)_LM analogously.
"""
import math, cmath
import mpmath as mp
from sympy.physics.wigner import wigner_9j, wigner_3j, clebsch_gordan
mp.mp.dps = 22
pi = mp.pi
def Y(l,m,O): return complex(mp.spherharm(l,m,O[0],O[1])) if abs(m)<=l else 0j
_cg={}
def CG(j1,m1,j2,m2,j3,m3):
    k=(j1,m1,j2,m2,j3,m3)
    if k not in _cg: _cg[k]=float(clebsch_gordan(j1,j2,j3,m1,m2,m3))
    return _cg[k]
_3j={}
def T3(l1,l2,l3,m1,m2,m3):
    k=(l1,l2,l3,m1,m2,m3)
    if k not in _3j: _3j[k]=float(wigner_3j(l1,l2,l3,m1,m2,m3))
    return _3j[k]
_9j={}
def N9(a,b,c,d,e,f,g,h,i):
    k=(a,b,c,d,e,f,g,h,i)
    if k not in _9j: _9j[k]=float(wigner_9j(a,b,c,d,e,f,g,h,i))
    return _9j[k]
def report(tag,w,tol=1e-11):
    w=float(w); okk=w<tol; print(f"  {tag:34s} {'PASS' if okk else 'FAIL'}  worst={w:.2e}"); return okk
ok=True

def Bip(l1,l2,L,M,O1,O2):
    s=0j
    for m1 in range(-l1,l1+1):
        m2=M-m1
        if abs(m2)>l2: continue
        s+=CG(l1,m1,l2,m2,L,M)*Y(l1,m1,O1)*Y(l2,m2,O2)
    return s
def Trip(l1,l2,l3,lam,L,M,O1,O2,O3):  # {Y_l1 x {Y_l2 x Y_l3}_lam}_LM
    s=0j
    for m1 in range(-l1,l1+1):
        m23=M-m1
        if abs(m23)>lam: continue
        for m2 in range(-l2,l2+1):
            m3=m23-m2
            if abs(m3)>l3: continue
            s+=CG(l1,m1,lam,m23,L,M)*CG(l2,m2,l3,m3,lam,m23)*Y(l1,m1,O1)*Y(l2,m2,O2)*Y(l3,m3,O3)
    return s

O1=(0.7,1.1); O2=(1.9,4.2); O3=(2.3,0.5)

print("Sec 5.16 bipolar")
# 5.16.4  sum_{LM} |Bip|^2 = (2l1+1)(2l2+1)/(4pi)^2
w=0.0
for l1 in range(0,4):
    for l2 in range(0,4):
        s=sum(abs(Bip(l1,l2,L,M,O1,O2))**2 for L in range(abs(l1-l2),l1+l2+1) for M in range(-L,L+1))
        w=max(w,abs(s-(2*l1+1)*(2*l2+1)/(4*pi)**2))
ok&=report("5.16.4 sum|Bip|^2", w)
# 5.16.9 + 5.16.10  Bip_00 = (-1)^l1/sqrt(2l1+1) sum_m Y*_{l1 m}(O1)Y_{l1 m}(O2) delta_{l1 l2}
w=0.0
for l1 in range(0,4):
    for l2 in range(0,4):
        lhs=Bip(l1,l2,0,0,O1,O2)
        sp=sum(Y(l1,m,O1).conjugate()*Y(l1,m,O2) for m in range(-l1,l1+1))
        rhs=((-1)**l1/math.sqrt(2*l1+1)*sp) if l1==l2 else 0j
        w=max(w,abs(lhs-rhs))
ok&=report("5.16.9/10 Bip_00 scalar product", w)
# 5.16.6  inversion: Bip at (pi-th,ph+pi) = (-1)^{l1+l2} Bip
def inv(O): return (math.pi-O[0], O[1]+math.pi)
w=0.0
for l1 in range(0,3):
    for l2 in range(0,3):
        for L in range(abs(l1-l2),l1+l2+1):
            for M in range(-L,L+1):
                w=max(w,abs(Bip(l1,l2,L,M,inv(O1),inv(O2))-(-1)**(l1+l2)*Bip(l1,l2,L,M,O1,O2)))
ok&=report("5.16.6 inversion phase", w)
# 5.16.7/8  product of two bipolars = sum_LM C^{LM}_{L'M'L''M''} sum_{l1l2} B^{l1l2L}_{...} Bip(l1,l2,L,M)
#  B = sqrt((2l1'+1)(2l2'+1)(2l1''+1)(2l2''+1)(2L'+1)(2L''+1)/(4pi)^2)
#      C^{l1 0}_{l1'0 l1''0} C^{l2 0}_{l2'0 l2''0} { l1' l1'' l1 ; l2' l2'' l2 ; L' L'' L }
def Bcoef(l1p,l2p,Lp,l1pp,l2pp,Lpp,l1,l2,L):
    return (math.sqrt((2*l1p+1)*(2*l2p+1)*(2*l1pp+1)*(2*l2pp+1)*(2*Lp+1)*(2*Lpp+1)/(4*pi)**2)
            *CG(l1p,0,l1pp,0,l1,0)*CG(l2p,0,l2pp,0,l2,0)
            *N9(l1p,l1pp,l1, l2p,l2pp,l2, Lp,Lpp,L))
w=0.0
for (l1p,l2p,l1pp,l2pp) in [(1,1,1,1),(2,1,1,2),(2,2,1,1)]:
    for Lp in range(abs(l1p-l2p),l1p+l2p+1):
        for Lpp in range(abs(l1pp-l2pp),l1pp+l2pp+1):
            for Mp in range(-Lp,Lp+1):
                for Mpp in range(-Lpp,Lpp+1):
                    lhs=Bip(l1p,l2p,Lp,Mp,O1,O2)*Bip(l1pp,l2pp,Lpp,Mpp,O1,O2)
                    s=0j
                    for L in range(abs(Lp-Lpp),Lp+Lpp+1):
                        M=Mp+Mpp
                        if abs(M)>L: continue
                        for l1 in range(abs(l1p-l1pp),l1p+l1pp+1):
                            for l2 in range(abs(l2p-l2pp),l2p+l2pp+1):
                                if L<abs(l1-l2) or L>l1+l2: continue
                                s+=CG(Lp,Mp,Lpp,Mpp,L,M)*Bcoef(l1p,l2p,Lp,l1pp,l2pp,Lpp,l1,l2,L)*Bip(l1,l2,L,M,O1,O2)
                    w=max(w,abs(lhs-s))
ok&=report("5.16.7/8 product (9j)", w)

print("\nSec 5.16 tripolar")
# 5.16.14  sum_{lam L M} |Trip|^2 = (2l1+1)(2l2+1)(2l3+1)/(4pi)^3
w=0.0
for l1 in range(0,3):
    for l2 in range(0,3):
        for l3 in range(0,3):
            s=0.0
            for lam in range(abs(l2-l3),l2+l3+1):
                for L in range(abs(l1-lam),l1+lam+1):
                    for M in range(-L,L+1):
                        s+=abs(Trip(l1,l2,l3,lam,L,M,O1,O2,O3))**2
            w=max(w,abs(s-(2*l1+1)*(2*l2+1)*(2*l3+1)/(4*pi)**3))
ok&=report("5.16.14 sum|Trip|^2", w)
# 5.16.19  Trip_{lam,00} = (-1)^{l1+l2+l3} delta_{lam,l1} sum 3j(l1l2l3;m1m2m3) Y Y Y
w=0.0
for l1 in range(0,3):
    for l2 in range(0,3):
        for l3 in range(abs(l1-l2),l1+l2+1):
            for lam in range(abs(l2-l3),l2+l3+1):
                lhs=Trip(l1,l2,l3,lam,0,0,O1,O2,O3)
                if lam==l1:
                    rhs=0j
                    for m1 in range(-l1,l1+1):
                        for m2 in range(-l2,l2+1):
                            m3=-(m1+m2)
                            if abs(m3)>l3: continue
                            rhs+=T3(l1,l2,l3,m1,m2,m3)*Y(l1,m1,O1)*Y(l2,m2,O2)*Y(l3,m3,O3)
                    rhs*=(-1)**(l1+l2+l3)
                else:
                    rhs=0j
                w=max(w,abs(lhs-rhs))
ok&=report("5.16.19 Trip_00 (3j)", w)

# 5.16.17/18  tripolar CG series (B coefficient carries TWO 9j symbols)
def Bt(l1p,l2p,l3p,lamp,Lp, l1pp,l2pp,l3pp,lampp,Lpp, l1,l2,l3,lam,L):
    return (math.sqrt((2*l1p+1)*(2*l1pp+1)*(2*l2p+1)*(2*l2pp+1)*(2*l3p+1)*(2*l3pp+1)*(2*Lp+1)*(2*Lpp+1)/(4*pi)**3)
            *math.sqrt((2*lamp+1)*(2*lampp+1)*(2*lam+1))
            *CG(l1p,0,l1pp,0,l1,0)*CG(l2p,0,l2pp,0,l2,0)*CG(l3p,0,l3pp,0,l3,0)
            *N9(l1p,l1pp,l1, lamp,lampp,lam, Lp,Lpp,L)
            *N9(l2p,l2pp,l2, l3p,l3pp,l3, lamp,lampp,lam))
w=0.0
(l1p,l2p,l3p)=(1,1,0); (l1pp,l2pp,l3pp)=(1,0,1)
for lamp in range(abs(l2p-l3p),l2p+l3p+1):
  for lampp in range(abs(l2pp-l3pp),l2pp+l3pp+1):
    for Lp in range(abs(l1p-lamp),l1p+lamp+1):
      for Lpp in range(abs(l1pp-lampp),l1pp+lampp+1):
        for Mp in range(-Lp,Lp+1):
          for Mpp in range(-Lpp,Lpp+1):
            lhs=Trip(l1p,l2p,l3p,lamp,Lp,Mp,O1,O2,O3)*Trip(l1pp,l2pp,l3pp,lampp,Lpp,Mpp,O1,O2,O3)
            s=0j
            for L in range(abs(Lp-Lpp),Lp+Lpp+1):
                M=Mp+Mpp
                if abs(M)>L: continue
                for l1 in range(abs(l1p-l1pp),l1p+l1pp+1):
                    for l2 in range(abs(l2p-l2pp),l2p+l2pp+1):
                        for l3 in range(abs(l3p-l3pp),l3p+l3pp+1):
                            for lam in range(abs(l2-l3),l2+l3+1):
                                if L<abs(l1-lam) or L>l1+lam: continue
                                s+=CG(Lp,Mp,Lpp,Mpp,L,M)*Bt(l1p,l2p,l3p,lamp,Lp,l1pp,l2pp,l3pp,lampp,Lpp,l1,l2,l3,lam,L)*Trip(l1,l2,l3,lam,L,M,O1,O2,O3)
            w=max(w,abs(lhs-s))
ok&=report("5.16.17/18 tripolar product (2x9j)", w)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
