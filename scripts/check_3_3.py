#!/usr/bin/env python3
r"""Verify Sec 3.3 recoupling relations for irreducible tensor products.

Irreducible tensor of rank a = array {M: complex}, M=-a..a.
Coupled product {P_a (x) Q_b}_c : [.]_M = sum_{M1M2} C^{cM}_{a M1,b M2} P[M1] Q[M2].
Scalar product (book 3.1.30): (X_J . Y_J) = sum_M (-1)^{-M} X_M Y_{-M}
   = (-1)^{-J} sqrt(2J+1) {X (x) Y}_00.
Pi factor: Pi_{ab..} = sqrt(prod (2x+1)).  6j via Racah; 9j via 6j sum.
"""
import cmath, random
random.seed(7)
from math import sqrt as msqrt
import mpmath as mp
mp.mp.dps=30

def _tri(a,b,c): return abs(a-b)<=c<=a+b and (a+b+c)==int(a+b+c)
def cg(j1,m1,j2,m2,j3,m3):
    from sympy.physics.quantum.cg import CG
    from sympy import Rational as R, N
    if m1+m2!=m3 or not _tri(j1,j2,j3): return 0.0
    return float(N(CG(R(j1),R(m1),R(j2),R(m2),R(j3),R(m3)).doit(),25))
def sixj(j1,j2,j3,j4,j5,j6):
    tris=[(j1,j2,j3),(j1,j5,j6),(j4,j2,j6),(j4,j5,j3)]
    if not all(_tri(*t) for t in tris): return 0.0
    def dl(a,b,c):
        from math import factorial as f
        return msqrt(f(int(a+b-c))*f(int(a-b+c))*f(int(-a+b+c))/f(int(a+b+c+1)))
    from math import factorial as f
    pref=dl(j1,j2,j3)*dl(j1,j5,j6)*dl(j4,j2,j6)*dl(j4,j5,j3)
    a=[j1+j2+j3,j1+j5+j6,j4+j2+j6,j4+j5+j3]; b=[j1+j2+j4+j5,j2+j3+j5+j6,j1+j3+j4+j6]
    s=0.0
    for t in range(int(max(a)),int(min(b))+1):
        s+=(-1)**t*f(t+1)/(f(int(t-a[0]))*f(int(t-a[1]))*f(int(t-a[2]))*f(int(t-a[3]))
             *f(int(b[0]-t))*f(int(b[1]-t))*f(int(b[2]-t)))
    return pref*s
def ninej(a,b,c,d,e,f,g,h,k):
    s=0.0
    xs=[x for x in range(int(max(abs(a-k),abs(b-f),abs(d-h))),int(min(a+k,b+f,d+h))+1)]
    for x in xs:
        s+=(-1)**(2*x)*(2*x+1)*sixj(a,b,c,f,k,x)*sixj(d,e,f,b,x,h)*sixj(g,h,k,x,a,d)
    return s
def Pi(*js): return msqrt(1.0)*mp.sqrt(mp.fprod([2*j+1 for j in js]))

def rt(a):
    return {M: complex(random.uniform(-1,1),random.uniform(-1,1)) for M in range(-a,a+1)}
def cpl(P,a,Q,b,c):
    out={}
    for M in range(-c,c+1):
        s=0.0
        for M1 in range(-a,a+1):
            M2=M-M1
            if -b<=M2<=b: s+=cg(a,M1,b,M2,c,M)*P[M1]*Q[M2]
        out[M]=s
    return out
def sc(X,Y,J):  # (X.Y) book 3.1.30
    return sum((-1)**(-M)*X[M]*Y[-M] for M in range(-J,J+1))
def dictnorm(A,B): return max(abs(A[M]-B[M]) for M in A)
def report(tag,w,tol=1e-9):
    w=float(abs(w)); ok=w<tol; print(f"  {tag:44s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok
ok=True

# random tensors
A=rt(1);B=rt(1);C=rt(1);D=rt(1)
Pa,Qb,Rd,Se={},{},{},{}

print("sanity: 3.1.34 (X.Y)=(-1)^-J sqrt(2J+1){XY}_00")
for J in (1,2):
    X=rt(J);Y=rt(J)
    t00=cpl(X,J,Y,J,0)[0]
    ok&=report(f"J={J} scalar=rank0", abs(sc(X,Y,J)-((-1)**(-J))*mp.sqrt(2*J+1)*t00))

print("\n3.3.1 {{PQ}c R}f = (-1)^{a+b+f+d} sum_h Pi_hc 6j {..} {P{QR}h}f")
def rel331(a,b,c,d,f):
    P=rt(a);Q=rt(b);R=rt(d)
    lhs=cpl(cpl(P,a,Q,b,c),c,R,d,f)
    rhs={M:0.0 for M in range(-f,f+1)}
    for h in range(abs(b-d),b+d+1):
        coef=(-1)**(a+b+f+d)*Pi(h,c)*sixj(a,b,c,d,f,h)
        term=cpl(P,a,cpl(Q,b,R,d,h),h,f)
        for M in rhs: rhs[M]+=coef*term[M]
    return dictnorm(lhs,rhs)
w=0.0
for (a,b,c,d,f) in [(1,1,2,1,1),(1,1,0,1,1),(2,1,2,1,2),(1,2,3,2,1),(1,1,1,1,2)]:
    w=max(w,rel331(a,b,c,d,f))
ok&=report("3.3.1", w)

print("3.3.7 {P{QR}f}e symmetry relations")
def rel337(a,b,d,f,e):
    P=rt(a);Q=rt(b);R=rt(d)
    lhs=cpl(P,a,cpl(Q,b,R,d,f),f,e)
    # (-1)^{b+d-f}{P{RQ}f}e
    r1=cpl(P,a,cpl(R,d,Q,b,f),f,e)
    v1=max(abs(lhs[M]-(-1)**(b+d-f)*r1[M]) for M in lhs)
    # (-1)^{a+f-e}{{QR}f P}e
    r2=cpl(cpl(Q,b,R,d,f),f,P,a,e)
    v2=max(abs(lhs[M]-(-1)**(a+f-e)*r2[M]) for M in lhs)
    return max(v1,v2)
w=0.0
for (a,b,d,f,e) in [(1,1,1,1,1),(1,1,1,2,1),(2,1,1,2,2),(1,2,1,1,1)]:
    w=max(w,rel337(a,b,d,f,e))
ok&=report("3.3.7", w)

print("3.3.8 {{PQ}c R}f = (-1)^{c+d+f} sum_h Pi_ch 6j {Q{PR}h}f  (commuting)")
def rel338(a,b,c,d,f):
    P=rt(a);Q=rt(b);R=rt(d)
    lhs=cpl(cpl(P,a,Q,b,c),c,R,d,f)
    rhs={M:0.0 for M in range(-f,f+1)}
    for h in range(abs(a-d),a+d+1):
        coef=(-1)**(c+d+f)*Pi(c,h)*sixj(a,b,c,f,d,h)
        term=cpl(Q,b,cpl(P,a,R,d,h),h,f)
        for M in rhs: rhs[M]+=coef*term[M]
    return dictnorm(lhs,rhs)
w=0.0
for (a,b,c,d,f) in [(1,1,2,1,1),(2,1,2,1,2),(1,2,2,1,1),(1,1,1,1,0)]:
    w=max(w,rel338(a,b,c,d,f))
ok&=report("3.3.8", w)

print("3.3.11 9j recoupling of 4 tensors")
def rel3311(a,b,c,d,e,f,k):
    P=rt(a);Q=rt(b);R=rt(d);S=rt(e)
    lhs=cpl(cpl(P,a,Q,b,c),c,cpl(R,d,S,e,f),f,k)
    rhs={M:0.0 for M in range(-k,k+1)}
    for g in range(abs(a-d),a+d+1):
        for h in range(abs(b-e),b+e+1):
            coef=Pi(c,f,g,h)*ninej(a,b,c,d,e,f,g,h,k)
            term=cpl(cpl(P,a,R,d,g),g,cpl(Q,b,S,e,h),h,k)
            for M in rhs: rhs[M]+=coef*term[M]
    return dictnorm(lhs,rhs)
w=0.0
for (a,b,c,d,e,f,k) in [(1,1,1,1,1,1,1),(1,1,2,1,1,2,1),(1,1,0,1,1,0,0),(1,1,1,1,1,2,2)]:
    w=max(w,rel3311(a,b,c,d,e,f,k))
ok&=report("3.3.11", w)

print("3.3.13 scalar ({PQ}c.{RS}c) = (-1)^{2a+b-d} sum_g Pi_c^2 6j (..)")
def rel3313(a,b,c,d,e):
    P=rt(a);Q=rt(b);R=rt(d);S=rt(e)
    L1=cpl(P,a,Q,b,c); L2=cpl(R,d,S,e,c)
    lhs=sc(L1,L2,c)
    rhs=0.0
    for g in range(max(abs(a-d),abs(b-e)),min(a+d,b+e)+1):
        rhs+=(-1)**(2*a+b-d)*Pi(c)**2*sixj(a,b,c,e,d,g)*sc(cpl(P,a,R,d,g),cpl(Q,b,S,e,g),g)
    return abs(lhs-rhs)
w=0.0
for (a,b,c,d,e) in [(1,1,1,1,1),(1,1,2,1,1),(2,1,2,1,2),(1,2,1,2,1)]:
    w=max(w,rel3313(a,b,c,d,e))
ok&=report("3.3.13", w)

print("3.3.2 ({PQ}c.R_c) = (-1)^{-c+a} Pi_c/Pi_a (P_a.{QR}_a)")
def rel332(a,b,c):
    P=rt(a);Q=rt(b);R=rt(c)
    lhs=sc(cpl(P,a,Q,b,c),R,c)
    rhs=(-1)**(-c+a)*Pi(c)/Pi(a)*sc(P,cpl(Q,b,R,c,a),a)
    return abs(lhs-rhs)
w=0.0
for (a,b,c) in [(1,1,1),(1,1,2),(2,1,2),(1,2,2),(2,2,1)]:
    w=max(w,rel332(a,b,c))
ok&=report("3.3.2", w)

print("3.3.9 ({PQ}c.R_c) = (-1)^{-a} Pi_c/Pi_b (Q_b.{PR}_b)  (commuting)")
def rel339(a,b,c):
    P=rt(a);Q=rt(b);R=rt(c)
    lhs=sc(cpl(P,a,Q,b,c),R,c)
    rhs=(-1)**(-a)*Pi(c)/Pi(b)*sc(Q,cpl(P,a,R,c,b),b)
    return abs(lhs-rhs)
w=0.0
for (a,b,c) in [(1,1,1),(1,1,2),(2,1,2),(1,2,2),(2,2,2)]:
    w=max(w,rel339(a,b,c))
ok&=report("3.3.9", w)

print("3.3.12 ({PQ}c {RS}f)k via 6j x 6j = {{{PR}h Q}g S}k")
def rel3312(a,b,c,d,e,f,k):
    P=rt(a);Q=rt(b);R=rt(d);S=rt(e)
    lhs=cpl(cpl(P,a,Q,b,c),c,cpl(R,d,S,e,f),f,k)
    rhs={M:0.0 for M in range(-k,k+1)}
    for g in range(abs(k-e),k+e+1):
        for h in range(abs(a-d),a+d+1):
            coef=(-1)**(h+b-k-e)*Pi(c,f,g,h)*sixj(a,b,c,g,d,h)*sixj(d,e,f,k,c,g)  # book Pi^2 is a misprint
            term=cpl(cpl(cpl(P,a,R,d,h),h,Q,b,g),g,S,e,k)
            for M in rhs: rhs[M]+=coef*term[M]
    return dictnorm(lhs,rhs)
w=0.0
for (a,b,c,d,e,f,k) in [(1,1,1,1,1,1,1),(1,1,2,1,1,1,1),(1,1,1,1,1,2,1)]:
    w=max(w,rel3312(a,b,c,d,e,f,k))
ok&=report("3.3.12", w)

print("3.3.14 ({PQ}c.{RS}c) = sum_h (-1)^{e+b+d+h} Pi_cch/Pi_e 6j ({{PR}h Q}e . S)")
def rel3314(a,b,c,d,e):
    P=rt(a);Q=rt(b);R=rt(d);S=rt(e)
    lhs=sc(cpl(P,a,Q,b,c),cpl(R,d,S,e,c),c)
    rhs=0.0
    for h in range(abs(a-d),a+d+1):
        rhs+=(-1)**(e+b+d+h)*Pi(c,c,h)/Pi(e)*sixj(a,b,c,e,d,h)*sc(cpl(cpl(P,a,R,d,h),h,Q,b,e),S,e)
    return abs(lhs-rhs)
w=0.0
for (a,b,c,d,e) in [(1,1,1,1,1),(1,1,2,1,1),(2,1,2,1,2)]:
    w=max(w,rel3314(a,b,c,d,e))
ok&=report("3.3.14", w)

# ---- non-commuting relations 3.3.16-3.3.23 (matrix reps) ----
print("\n3.3.16-3.3.23 non-commuting tensors (matrix reps)")
import numpy as np
def rtM(a,dim=3):
    return {M: (np.random.randn(dim,dim)+1j*np.random.randn(dim,dim)) for M in range(-a,a+1)}
def cplM(P,a,Q,b,c):
    out={}
    for M in range(-c,c+1):
        s=np.zeros_like(next(iter(P.values())))
        for M1 in range(-a,a+1):
            M2=M-M1
            if -b<=M2<=b: s=s+cg(a,M1,b,M2,c,M)*P[M1]@Q[M2]
        out[M]=s
    return out
def commutM(P,a,Q,b,c):  # R^{ab}_c = sum C^{cM}_{aM1,bM2}[P_M1,Q_M2]
    out={}
    for M in range(-c,c+1):
        s=np.zeros_like(next(iter(P.values())))
        for M1 in range(-a,a+1):
            M2=M-M1
            if -b<=M2<=b: s=s+cg(a,M1,b,M2,c,M)*(P[M1]@Q[M2]-Q[M2]@P[M1])
        out[M]=s
    return out
def dnormM(A,B): return max(np.max(np.abs(A[M]-B[M])) for M in A)
np.random.seed(3)
# 3.3.16: {{PQ}c R}f = (-1)^{f+d+c} sum_h Pi_ch 6j {Q{PR}h}f + {R^{ab}_c R_d}f
def rel3316(a,b,c,d,f):
    P=rtM(a);Q=rtM(b);R=rtM(d)
    lhs=cplM(cplM(P,a,Q,b,c),c,R,d,f)
    rhs={M:np.zeros_like(next(iter(P.values()))) for M in range(-f,f+1)}
    for h in range(abs(a-d),a+d+1):
        coef=(-1)**(f+d+c)*Pi(c,h)*sixj(a,b,c,f,d,h)
        term=cplM(Q,b,cplM(P,a,R,d,h),h,f)
        for M in rhs: rhs[M]=rhs[M]+float(coef)*term[M]
    Rc=commutM(P,a,Q,b,c)
    term2=cplM(Rc,c,R,d,f)
    for M in rhs: rhs[M]=rhs[M]+term2[M]
    return dnormM(lhs,rhs)
w=0.0
for (a,b,c,d,f) in [(1,1,2,1,1),(1,1,1,1,1),(2,1,2,1,2)]:
    w=max(w,rel3316(a,b,c,d,f))
ok&=report("3.3.16 (non-commuting)", w, tol=1e-8)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
