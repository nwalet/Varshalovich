#!/usr/bin/env python3
r"""Verify Sec 3.2.1 vector<->irreducible-tensor identities (3.2.2-3.2.23).

Spherical comps of a vector: V_{+1}=-(Vx+iVy)/sqrt2, V0=Vz, V_{-1}=(Vx-iVy)/sqrt2.
{A_1 (x) B_1}_{lM} = sum_{mu nu} C^{lM}_{1 mu,1 nu} A_mu B_nu.
Coupled products of rank-1 tensors, compared to cartesian vector RHS.
"""
import numpy as np
from math import sqrt
np.random.seed(11)
def cg(j1,m1,j2,m2,j3,m3):
    from sympy.physics.quantum.cg import CG
    from sympy import Rational as R, N
    if m1+m2!=m3: return 0.0
    return float(N(CG(R(j1),R(m1),R(j2),R(m2),R(j3),R(m3)).doit(),25))
r2=sqrt(2)
def sph(V): return {1:-(V[0]+1j*V[1])/r2, 0:V[2]+0j, -1:(V[0]-1j*V[1])/r2}
def cart(comp):  # spherical comps dict -> cartesian 3-vector
    Vx=(comp[-1]-comp[1])/r2; Vy=1j*(comp[1]+comp[-1])/r2; Vz=comp[0]
    return np.array([Vx,Vy,Vz])
def coupleT(P,a,Q,b,c):  # both dicts of spherical comps
    out={}
    for M in range(-c,c+1):
        s=0j
        for m1 in range(-a,a+1):
            m2=M-m1
            if -b<=m2<=b: s+=cg(a,m1,b,m2,c,M)*P[m1]*Q[m2]
        out[M]=s
    return out
def dot(A,B): return complex(np.dot(A,B))
def cross(A,B): return np.cross(A,B)
def report(tag,w,tol=1e-11):
    w=float(abs(w)); ok=w<tol; print(f"  {tag:36s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok
def vdiff(comp,vec): return float(np.max(np.abs(cart(comp)-vec)))
ok=True
A=np.random.randn(3)+1j*np.random.randn(3); B=np.random.randn(3)+1j*np.random.randn(3)
C=np.random.randn(3)+1j*np.random.randn(3); D=np.random.randn(3)+1j*np.random.randn(3)
sA,sB,sC,sD=sph(A),sph(B),sph(C),sph(D)

print("3.2.2/4/6 basic products")
ok&=report("3.2.2 {AB}_00=-1/sqrt3 (A.B)", abs(coupleT(sA,1,sB,1,0)[0]-(-1/sqrt(3)*dot(A,B))))
c1=coupleT(sA,1,sB,1,1); ok&=report("3.2.4 {AB}_1=i/sqrt2 [AxB]", vdiff(c1,1j/r2*cross(A,B)))
# 3.2.6/7 explicit rank-2
c2=coupleT(sA,1,sB,1,2)
e22=sA[1]*sB[1]; e21=1/r2*(sA[1]*sB[0]+sA[0]*sB[1]); e20=1/sqrt(6)*(sA[1]*sB[-1]+2*sA[0]*sB[0]+sA[-1]*sB[1])
ok&=report("3.2.7 {AB}_2M explicit", max(abs(c2[2]-e22),abs(c2[1]-e21),abs(c2[0]-e20)))

print("3.2.8-11 triple products")
# 3.2.8 {{AB}_0 C}_1 = -1/sqrt3 (A.B) C
t=coupleT(coupleT(sA,1,sB,1,0),0,sC,1,1); ok&=report("3.2.8", vdiff(t,-1/sqrt(3)*dot(A,B)*C))
# 3.2.9 {{AB}_1 C}_0 = -i/sqrt6 [AxB].C
t=coupleT(coupleT(sA,1,sB,1,1),1,sC,1,0)[0]; ok&=report("3.2.9", abs(t-(-1j/sqrt(6)*dot(cross(A,B),C))))
# 3.2.10 {{AB}_1 C}_1 = 1/2 A(B.C) - 1/2 B(A.C)
t=coupleT(coupleT(sA,1,sB,1,1),1,sC,1,1); ok&=report("3.2.10", vdiff(t,0.5*A*dot(B,C)-0.5*B*dot(A,C)))
# 3.2.11 {{AB}_2 C}_1 = sqrt(3/5){1/3 C(A.B) - 1/2 B(A.C) - 1/2 A(B.C)}
t=coupleT(coupleT(sA,1,sB,1,2),2,sC,1,1)
ok&=report("3.2.11", vdiff(t,sqrt(3/5)*(1/3*C*dot(A,B)-0.5*B*dot(A,C)-0.5*A*dot(B,C))))

print("3.2.12-20 quadruple products")
def q(c1,c2):  # {{AB}_c1 {CD}_c2}
    return coupleT(coupleT(sA,1,sB,1,c1),c1,coupleT(sC,1,sD,1,c2),c2, None)
def qk(c1,c2,k):
    return coupleT(coupleT(sA,1,sB,1,c1),c1,coupleT(sC,1,sD,1,c2),c2,k)
# 3.2.12 {{AB}0{CD}0}0 = 1/3 (A.B)(D.C)
ok&=report("3.2.12", abs(qk(0,0,0)[0]-(1/3*dot(A,B)*dot(D,C))))
# 3.2.13 {{AB}1{CD}0}1 = -i/sqrt6 [AxB](C.D)
ok&=report("3.2.13", vdiff(qk(1,0,1),-1j/sqrt(6)*cross(A,B)*dot(C,D)))
# 3.2.14 {{AB}0{CD}1}1 = -i/sqrt6 (A.B)[CxD]
ok&=report("3.2.14", vdiff(qk(0,1,1),-1j/sqrt(6)*dot(A,B)*cross(C,D)))
# 3.2.15 {{AB}1{CD}1}0 = 1/(2sqrt3){(A.C)(B.D)-(A.D)(B.C)}
ok&=report("3.2.15", abs(qk(1,1,0)[0]-(1/(2*sqrt(3))*(dot(A,C)*dot(B,D)-dot(A,D)*dot(B,C)))))
# 3.2.16 {{AB}1{CD}1}1 = -i/(2sqrt2){C(D.[AxB])-D(C.[AxB])}
ok&=report("3.2.16", vdiff(qk(1,1,1),-1j/(2*r2)*(C*dot(D,cross(A,B))-D*dot(C,cross(A,B)))))
# 3.2.19 {{AB}2{CD}2}0 = 1/sqrt5{1/2(A.C)(B.D)-1/3(A.B)(C.D)+1/2(A.D)(B.C)}
ok&=report("3.2.19", abs(qk(2,2,0)[0]-(1/sqrt(5)*(0.5*dot(A,C)*dot(B,D)-1/3*dot(A,B)*dot(C,D)+0.5*dot(A,D)*dot(B,C)))))
# 3.2.20 {{AB}2{CD}2}1 = -i/(2sqrt10){(A.C)[BxD]+(A.D)[BxC]+(B.C)[AxD]+(B.D)[AxC]}
ok&=report("3.2.20", vdiff(qk(2,2,1),-1j/(2*sqrt(10))*(dot(A,C)*cross(B,D)+dot(A,D)*cross(B,C)+dot(B,C)*cross(A,D)+dot(B,D)*cross(A,C))))
# 3.2.17 {{AB}2{CD}1}1 = i sqrt3/sqrt10 {1/3(A.B)[CxD]-1/2 B(D.[AxC])-1/2 A(D.[BxC])}
ok&=report("3.2.17", vdiff(qk(2,1,1),1j*sqrt(3)/sqrt(10)*(1/3*dot(A,B)*cross(C,D)-0.5*B*dot(D,cross(A,C))-0.5*A*dot(D,cross(B,C)))))
# 3.2.18 {{AB}1{CD}2}1 = i sqrt3/sqrt10 {1/3(C.D)[AxB]-1/2 C(B.[DxA])-1/2 D(B.[CxA])}
ok&=report("3.2.18", vdiff(qk(1,2,1),1j*sqrt(3)/sqrt(10)*(1/3*dot(C,D)*cross(A,B)-0.5*C*dot(B,cross(D,A))-0.5*D*dot(B,cross(C,A)))))

print("3.2.23 {{...{AA}2 A}3...}_nm = sqrt(4pi n!/(2n+1)!!) |A|^n Y_nm")
import mpmath as mp
def Y(l,m,th,ph): return complex(mp.spherharm(l,m,th,ph))
# use a real unit-ish vector A0 with polar angles
th,ph=0.7,1.1; Amag=1.3
A0=Amag*np.array([mp.sin(th)*mp.cos(ph),mp.sin(th)*mp.sin(ph),mp.cos(th)],dtype=complex)
sA0=sph(A0)
def chain(n):
    # {{...{A A}_2 A}_3 ... A}_n
    cur=coupleT(sA0,1,sA0,1,2); rank=2
    for l in range(3,n+1):
        cur=coupleT(cur,rank,sA0,1,l); rank=l
    return cur
w=0.0
from math import factorial
def dfact(k):
    r=1
    while k>0: r*=k; k-=2
    return r
for n in (2,3,4):
    cur=chain(n)
    coef=sqrt(4*mp.pi*factorial(n)/dfact(2*n+1))*Amag**n
    for m in range(-n,n+1):
        w=max(w,abs(cur[m]-coef*Y(n,m,th,ph)))
ok&=report("3.2.23 n=2,3,4", w, tol=1e-10)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
