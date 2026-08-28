#!/usr/bin/env python3
r"""Verify Sec 5.17 scalar-function expansions (5.17.9-5.17.33) vs
mpmath.spherharm.  Key reduction: (Y_l(O1).Y_l(O2)) = (2l+1)/(4pi) P_l(cos w12)
[eq 5.17.9], so every scalar expansion becomes a 1-D Legendre identity in
x=cos w12 (and the radial variables).  5.17.9 itself is checked against the
true sum_m Y*_{lm}(O1) Y_{lm}(O2).
"""
import math, cmath
import mpmath as mp
mp.mp.dps = 30
pi = mp.pi
def Y(l,m,O): return complex(mp.spherharm(l,m,O[0],O[1])) if abs(m)<=l else 0j
def P(l,x): return mp.legendre(l,x)
def jl(l,x): return mp.sqrt(pi/(2*x))*mp.besselj(l+mp.mpf(1)/2,x)
def hl1(l,x): return jl(l,x)+1j*(mp.sqrt(pi/(2*x))*mp.bessely(l+mp.mpf(1)/2,x))
def dfac(n):
    n=int(round(n))
    if n<=0: return mp.mpf(1)
    r=mp.mpf(1)
    while n>1: r*=n; n-=2
    return r
def fac(n): return mp.factorial(int(round(n)))
def rf(a,l): return mp.rf(a,l)          # Pochhammer (a)_l
def F(a,b,c,z): return mp.hyp2f1(a,b,c,z)
def report(tag,w,tol=mp.mpf('1e-10')):
    w=float(w); ok=w<tol; print(f"  {tag:40s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok
ok=True
LMAX=200

# 5.17.9  (Y_l.Y_l) = sum_m Y*_{lm}(O1)Y_{lm}(O2) = (2l+1)/4pi P_l(cos w12)
def cosw(O1,O2):
    return math.cos(O1[0])*math.cos(O2[0])+math.sin(O1[0])*math.sin(O2[0])*math.cos(O1[1]-O2[1])
w=0.0
for O1 in [(0.7,1.1),(1.9,4.2),(2.3,0.5)]:
    for O2 in [(1.2,2.0),(0.4,5.1)]:
        x=cosw(O1,O2)
        for l in range(0,6):
            sp=sum(Y(l,m,O1).conjugate()*Y(l,m,O2) for m in range(-l,l+1))
            w=max(w, abs(sp-(2*l+1)/(4*pi)*P(l,x)))
ok&=report("5.17.9 addition theorem", w)

# helper: scalar product via P_l (validated above)
def SP(l,x): return (2*l+1)/(4*pi)*P(l,x)

XS=[mp.mpf('-0.7'),mp.mpf('-0.2'),mp.mpf('0.3'),mp.mpf('0.85')]
print("Sec 5.17.3  (r1.r2)^n")
# 5.17.10/11/12  n=1,2,3 explicit
def dot(r1,r2,x): return r1*r2*x
w=0.0
for r1,r2 in [(1.3,2.1)]:
    for x in XS:
        # 5.17.10
        w=max(w, abs(dot(r1,r2,x)-4*pi/3*r1*r2*SP(1,x)))
        # 5.17.11
        rhs=4*pi/3*r1**2*r2**2*(SP(0,x)+mp.mpf(2)/5*SP(2,x))
        w=max(w, abs(dot(r1,r2,x)**2-rhs))
        # 5.17.12
        rhs=4*pi/5*r1**3*r2**3*(SP(1,x)+mp.mpf(2)/7*SP(3,x))
        w=max(w, abs(dot(r1,r2,x)**3-rhs))
ok&=report("5.17.10/11/12 (r1.r2)^{1,2,3}", w)
# 5.17.13  general n:  (r1.r2)^n = 4pi r1^n r2^n sum_{l} n!/((n-l)!!(n+l+1)!!) SP(l,x)
w=0.0
for n in range(0,7):
    for r1,r2 in [(1.3,2.1)]:
        for x in XS:
            s=mp.mpf(0)
            for l in range(n%2, n+1, 2):
                s+=fac(n)/(dfac(n-l)*dfac(n+l+1))*SP(l,x)
            w=max(w, abs(dot(r1,r2,x)**n-4*pi*r1**n*r2**n*s))
ok&=report("5.17.13 (r1.r2)^n general", w)
# 5.17.14  e^{i(r1.r2)} = 4pi sum i^l j_l(r1 r2) SP(l,x)
w=0.0
for r1,r2 in [(0.8,1.1),(1.5,1.3)]:
    z=r1*r2
    for x in XS:
        s=sum((1j)**l*jl(l,z)*SP(l,x) for l in range(0,LMAX))
        w=max(w, abs(cmath.exp(1j*float(z*x))-4*pi*s))
ok&=report("5.17.14 exp(i r1.r2)", w)
# 5.17.16/17  power series f=sum c_n (r1.r2)^n -> f_l = 4pi sum_{n=l,l+2} c_n n!(r1r2)^n/((n-l)!!(n+l+1)!!)
# test with f(t)=cosh(t) => c_n = 1/n! for even n; check sum_l f_l SP(l,x) = cosh(r1 r2 x)
w=0.0
for r1,r2 in [(0.7,1.1)]:
    z=r1*r2
    for x in XS:
        tot=mp.mpf(0)
        for l in range(0,40):
            fl=mp.mpf(0)
            for n in range(l,60,2):
                cn=(mp.mpf(1)/fac(n) if n%2==0 else mp.mpf(0))   # cosh coefficients
                fl+=cn*fac(n)*z**n/(dfac(n-l)*dfac(n+l+1))
            fl*=4*pi
            tot+=fl*SP(l,x)
        w=max(w, abs(tot-mp.cosh(z*x)))
ok&=report("5.17.16/17 power-series f (cosh)", w)

print("Sec 5.17.4  functions of r=|r1-r2|")
def rmag(r1,r2,x): return mp.sqrt(r1*r1+r2*r2-2*r1*r2*x)
# 5.17.19  Helmholtz Green: e^{ikr}/r = 4pi i k sum j_l(k r1) h^1_l(k r2) SP(l,x), r1<r2
w=0.0
k=mp.mpf('1.3'); r1,r2=mp.mpf('0.7'),mp.mpf('1.9')
for x in XS:
    r=rmag(r1,r2,x)
    s=sum(jl(l,k*r1)*hl1(l,k*r2)*SP(l,x) for l in range(0,LMAX))
    w=max(w, abs(mp.e**(1j*k*r)/r-4*pi*1j*k*s))
ok&=report("5.17.19 Helmholtz Green (r1<r2)", w)
# 5.17.21  1/r = 4pi/r2 sum 1/(2l+1) (r1/r2)^l SP(l,x), r1<r2
w=0.0
for x in XS:
    r=rmag(r1,r2,x)
    s=sum(1/(2*l+1)*(r1/r2)**l*SP(l,x) for l in range(0,LMAX))
    w=max(w, abs(1/r-4*pi/r2*s))
ok&=report("5.17.21 1/r multipole (r1<r2)", w)
# 5.17.22  symmetric 1/r
w=0.0
for (r1,r2) in [(mp.mpf('0.7'),mp.mpf('1.9')),(mp.mpf('2.1'),mp.mpf('0.6'))]:
    D=r1*r1+r2*r2
    for x in XS:
        r=rmag(r1,r2,x)
        s=mp.mpf(0)
        for l in range(0,60):
            inner=sum(dfac(2*n-1)/(dfac(n-l)*dfac(l+n+1))*(r1*r2/D)**n for n in range(l,120,2))
            s+=inner*SP(l,x)
        w=max(w, abs(1/r-4*pi/mp.sqrt(D)*s))
ok&=report("5.17.22 1/r symmetric", w)
# 5.17.23/24/25  r, 1/r^3, 1/r^5  (r1<r2)
r1,r2=mp.mpf('0.7'),mp.mpf('1.9')
w=0.0
for x in XS:
    r=rmag(r1,r2,x)
    s=sum(1/(2*l+1)*(r1**l/r2**(l+1))*(r1**2/(2*l+3)-r2**2/(2*l-1))*SP(l,x) for l in range(0,LMAX))
    w=max(w, abs(r-4*pi*s))
ok&=report("5.17.23 r (r1<r2)", w)
w=0.0
for x in XS:
    r=rmag(r1,r2,x)
    s=sum((r1**l/r2**(l+1))*SP(l,x) for l in range(0,LMAX))
    w=max(w, abs(1/r**3-4*pi/(r2*r2-r1*r1)*s))
ok&=report("5.17.24 1/r^3 (r1<r2)", w)
w=0.0
for x in XS:
    r=rmag(r1,r2,x)
    s=sum((2*l-1)*(2*l+3)*(r1**l/r2**(l+1))*(r2**2/(2*l-1)-r1**2/(2*l+3))*SP(l,x) for l in range(0,LMAX))
    w=max(w, abs(1/r**5-4*pi/(3*(r2*r2-r1*r1)**3)*s))
ok&=report("5.17.25 1/r^5 (r1<r2)", w)

print("Sec 5.17.5  r^n via a_l^n hypergeometric")
# 5.17.26  r^n = 4pi sum 1/(2l+1) a_l^n SP(l,x) = sum a_l^n P_l(x)  (since SP=(2l+1)/4pi P_l)
# 5.17.27  a_l^n = (-n/2)_l/(1/2)_l r2^n (r1/r2)^l F(l-n/2, -1/2-n/2; l+3/2; r1^2/r2^2), r1<r2
def a27(l,n,r1,r2):
    return (rf(-mp.mpf(n)/2,l)/rf(mp.mpf(1)/2,l)*r2**n*(r1/r2)**l
            *F(l-mp.mpf(n)/2, -mp.mpf(1)/2-mp.mpf(n)/2, l+mp.mpf(3)/2, r1**2/r2**2))
w=0.0
for n in [1,2,3,-1,-3,-5,4]:
    for x in XS:
        r=rmag(r1,r2,x)
        s=sum(a27(l,n,r1,r2)*P(l,x) for l in range(0,LMAX))
        w=max(w, abs(r**n-s))
ok&=report("5.17.26/27 r^n via a_l^n [27]", w)
# 5.17.28  alternate a_l^n (r1<r2)
def a28(l,n,r1,r2):
    return (rf(-mp.mpf(n)/2,l)/rf(mp.mpf(1)/2,l)*r1**l*(r2*r2-r1*r1)**(n+2)/r2**(l+n+4)
            *F(l+2+mp.mpf(n)/2, mp.mpf(3)/2+mp.mpf(n)/2, l+mp.mpf(3)/2, r1**2/r2**2))
w=0.0
for n in [-1,-3,-5]:   # requires n+2 sensible; test negative odd
    for x in XS:
        r=rmag(r1,r2,x)
        s=sum(a28(l,n,r1,r2)*P(l,x) for l in range(0,LMAX))
        w=max(w, abs(r**n-s))
ok&=report("5.17.28 r^n via a_l^n [28]", w)
# 5.17.31/32  symmetric a_l^n
def a31(l,n,r1,r2):
    D=r1*r1+r2*r2
    return (rf(-mp.mpf(n)/2,l)/rf(mp.mpf(1)/2,l)*(r1*r2)**l/D**(l-mp.mpf(n)/2)
            *F(mp.mpf(l)/2-mp.mpf(n)/4, mp.mpf(l)/2-mp.mpf(n)/4+mp.mpf(1)/2, l+mp.mpf(3)/2, (2*r1*r2/D)**2))
def a32(l,n,r1,r2):
    return (rf(-mp.mpf(n)/2,l)/rf(mp.mpf(1)/2,l)*(r1*r2)**l/(r1+r2)**(2*l-n)
            *F(l-mp.mpf(n)/2, 1+l, 2*l+2, 4*r1*r2/(r1+r2)**2))
for name,af in (("5.17.31",a31),("5.17.32",a32)):
    w=0.0
    for (r1,r2) in [(mp.mpf('0.7'),mp.mpf('1.9')),(mp.mpf('2.1'),mp.mpf('0.6'))]:
        for n in [1,-1,-3,2]:
            for x in XS:
                r=rmag(r1,r2,x)
                s=sum(af(l,n,r1,r2)*P(l,x) for l in range(0,LMAX))
                w=max(w, abs(r**n-s))
    ok&=report(f"{name} symmetric a_l^n", w)
# 5.17.33  a_l^n = (-n/2)_l/(n/2+2)_l (r2^2-r1^2)^{n+2} a_l^{-n-4}
r1,r2=mp.mpf('0.7'),mp.mpf('1.9')
w=0.0
for n in [-1,-3,1]:
    for l in range(0,6):
        lhs=a27(l,n,r1,r2)
        rhs=rf(-mp.mpf(n)/2,l)/rf(mp.mpf(n)/2+2,l)*(r2*r2-r1*r1)**(n+2)*a27(l,-n-4,r1,r2)
        w=max(w, abs(lhs-rhs))
ok&=report("5.17.33 a_l^n / a_l^{-n-4} relation", w)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
