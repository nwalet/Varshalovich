#!/usr/bin/env python3
r"""Verify Sec 6.3 (spin functions for S=1) and emit correct action tables.

Spin-1 operators in the SPHERICAL basis (rows/cols ordered m=+1,0,-1), so that
chi_{1m} are the unit column vectors e_m:
  S_z=diag(1,0,-1); S_+ ladder; S_x,S_y from S_+/S_-; S_0=S_z,
  S_{+1}=-(S_x+iS_y)/sqrt2, S_{-1}=(S_x-iS_y)/sqrt2.
  T_{2M} = sum_{mu,nu} C^{2M}_{1mu,1nu} S_mu S_nu           (Eq 2.6.4)
  Q_{ik} from inverting 2.6.7 (T<->Q) + traceless condition.
Cartesian basis vectors chi_i (spherical rep, Eq 6.3.8):
  chi_x=(-1,0,1)/sqrt2, chi_y=i(1,0,1)/sqrt2, chi_z=(0,1,0).
clebsch{a}{b}{c}{d}{e}{f}=C^{e f}_{a b,c d}.
"""
import mpmath as mp
mp.mp.dps=30
I3=mp.eye(3)
# spherical-basis spin-1 matrices, order |+1>,|0>,|-1>
Sz=mp.matrix([[1,0,0],[0,0,0],[0,0,-1]])
r2=mp.sqrt(2)
Sp=mp.matrix([[0,r2,0],[0,0,r2],[0,0,0]])   # S_+ |0>=sqrt2|+1>, |-1>=sqrt2|0>
Sm=Sp.T
Sx=(Sp+Sm)/2; Sy=(Sp-Sm)/(2j)
Sph={0:Sz, 1:-(Sx+1j*Sy)/r2, -1:(Sx-1j*Sy)/r2}
Scart={'x':Sx,'y':Sy,'z':Sz}
def cg(a,b,c,d,e,f):
    from sympy.physics.quantum.cg import CG
    from sympy import Rational as R, N
    from fractions import Fraction
    def rr(x): return R(Fraction(float(x)).limit_denominator(8))
    return complex(N(CG(rr(a),rr(b),rr(c),rr(d),rr(e),rr(f)).doit()))
# T_2M via 2.6.4
Tq={}
for M in range(-2,3):
    acc=mp.zeros(3,3)
    for mu in (1,0,-1):
        for nu in (1,0,-1):
            acc+=cg(1,mu,1,nu,2,M)*(Sph[mu]*Sph[nu])
    Tq[M]=acc
# Q_ik from T (invert 2.6.7) + traceless
Qzz=mp.sqrt(mp.mpf(2)/3)*Tq[0]
QxxmQyy=Tq[2]+Tq[-2]; Qxy=(Tq[2]-Tq[-2])/(2j)
Qxz=(Tq[-1]-Tq[1])/2; Qyz=(-Tq[1]-Tq[-1])/(2j)
Qxx=(-Qzz+QxxmQyy)/2; Qyy=(-Qzz-QxxmQyy)/2
Q={'xx':Qxx,'yy':Qyy,'zz':Qzz,'xy':Qxy,'yx':Qxy,'xz':Qxz,'zx':Qxz,'yz':Qyz,'zy':Qyz}
# basis vectors
e={1:mp.matrix([1,0,0]),0:mp.matrix([0,1,0]),-1:mp.matrix([0,0,1])}   # chi_{1m}
chi={'x':mp.matrix([-1,0,1])/r2,'y':1j*mp.matrix([1,0,1])/r2,'z':mp.matrix([0,1,0])}
def mnorm(M): return max(abs(M[i]) for i in range(len(M))) if min(M.rows,M.cols)==1 else max(abs(M[i,j]) for i in range(M.rows) for j in range(M.cols))
def report(tag,w,tol=mp.mpf('1e-12')):
    w=float(w); ok=w<tol; print(f"  {tag:44s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok
ok=True

# ---- anchor: 6.3.1 eigen ----
print("6.3.1 eigenvalues")
S2=Sx*Sx+Sy*Sy+Sz*Sz
ok&=report("S^2=2 I", mnorm(S2-2*I3))
w=max(mnorm(Sz*e[m]-m*e[m]) for m in (1,0,-1)); ok&=report("S_z chi_1m=m chi", w)

# ---- 6.3.2/8 cartesian<->spherical basis vectors ----
print("\n6.3.2 relations chi_1m <-> chi_i")
ok&=report("chi_11=-(chi_x+i chi_y)/sqrt2", mnorm(e[1]-(-(chi['x']+1j*chi['y'])/r2)))
ok&=report("chi_1-1=(chi_x-i chi_y)/sqrt2", mnorm(e[-1]-((chi['x']-1j*chi['y'])/r2)))
ok&=report("chi_10=chi_z", mnorm(e[0]-chi['z']))

# ---- 6.3.20/21 product expansions (definitional anchor for T,Q normalization) ----
print("\n6.3.20/21 product expansions")
w=mp.mpf(0)
for m in (1,0,-1):
    for mpp in (1,0,-1):
        lhs=e[m]*e[mpp].H
        rhs=mp.mpf(1)/3*(1 if m==mpp else 0)*I3
        for mu in (1,0,-1): rhs=rhs+1/r2*cg(1,mpp,1,mu,1,m)*Sph[mu]
        for M in range(-2,3): rhs=rhs+mp.sqrt(mp.mpf(5)/3)*cg(1,mpp,2,M,1,m)*Tq[M]
        w=max(w,mnorm(lhs-rhs))
ok&=report("6.3.20 chi_1m chi_1m'^dag", w, tol=mp.mpf('1e-13'))
eps={('x','y','z'):1,('y','z','x'):1,('z','x','y'):1,('x','z','y'):-1,('z','y','x'):-1,('y','x','z'):-1}
def eps_ikl(i,k,l): return eps.get((i,k,l),0)
w=mp.mpf(0)
for i in 'xyz':
    for k in 'xyz':
        lhs=chi[i]*chi[k].H
        rhs=mp.mpf(1)/3*(1 if i==k else 0)*I3
        for l in 'xyz':
            if eps_ikl(i,k,l): rhs=rhs+1j/2*eps_ikl(i,k,l)*Scart[l]
        rhs=rhs-Q[i+k]
        w=max(w,mnorm(lhs-rhs))
ok&=report("6.3.21 chi_i chi_k^dag", w, tol=mp.mpf('1e-13'))

# ---- 6.3.24 spin action (anchor) ----
print("\n6.3.24 S_mu chi_1m, S_i chi_k")
w=mp.mpf(0)
for mu in (1,0,-1):
    for m in (1,0,-1):
        lhs=Sph[mu]*e[m]; rhs=mp.zeros(3,1)
        for mpp in (1,0,-1): rhs=rhs+r2*cg(1,m,1,mu,1,mpp)*e[mpp]
        w=max(w,mnorm(lhs-rhs))
ok&=report("6.3.24 S_mu chi_1m", w, tol=mp.mpf('1e-13'))
w=mp.mpf(0)
for i in 'xyz':
    for k in 'xyz':
        lhs=Scart[i]*chi[k]; rhs=mp.zeros(3,1)
        for l in 'xyz':
            if eps_ikl(i,k,l): rhs=rhs+1j*eps_ikl(i,k,l)*chi[l]
        w=max(w,mnorm(lhs-rhs))
ok&=report("6.3.24 S_i chi_k = i eps chi_l", w)

# ---- emit correct action tables for the garbled detail eqs ----
def fmt(v):
    # identify v as combination of basis vectors, print compactly
    return v
def emit_actions(name, ops, basis, baslabel):
    print(f"\n== correct {name} ==")
    for oi,op in ops.items():
        row=[]
        for bl,bv in basis.items():
            res=op*bv
            row.append(f"{oi}.{bl}="+vecstr(res,basis))
        print("  "+" | ".join(row))
def vecstr(res,basis):
    # express res in terms of basis vectors
    terms=[]
    for bl,bv in basis.items():
        c=(bv.H*res)[0]/(bv.H*bv)[0]
        if abs(c)>1e-12:
            terms.append(f"({mp.nstr(mp.chop(c),4)}){bl}")
    return "+".join(terms) if terms else "0"

emit_actions("6.3.25 S_mu chi_1m", Sph, e, 'm')
emit_actions("6.3.26 S_mu chi_i", Sph, chi, 'i')
emit_actions("6.3.27 S_i chi_1m", Scart, e, 'm')
emit_actions("6.3.28 S_i chi_i", Scart, chi, 'i')
emit_actions("6.3.30 T_2M chi_1m", Tq, e, 'm')
emit_actions("6.3.31 T_2M chi_i", Tq, chi, 'i')
Qops={k:Q[k] for k in ('xx','yy','zz','xy','xz','yz')}
emit_actions("6.3.32 Q_ik chi_1m", Qops, e, 'm')
emit_actions("6.3.33 Q_ik chi_i", Qops, chi, 'i')

# ---- 6.3.23 emit cartesian decompositions chi_i chi_k^dag = 1/3 d I + i/2 eps S - Q ----
print("\n== correct 6.3.23 chi_i chi_k^dag (cartesian RHS) ==")
for i in 'xyz':
    for k in 'xyz':
        terms=[]
        if i==k: terms.append("I/3")
        for l in 'xyz':
            if eps_ikl(i,k,l): terms.append(f"{'+' if eps_ikl(i,k,l)>0 else '-'}i/2 S_{l}")
        terms.append(f"-Q_{i}{k}")
        print(f"  chi_{i} chi_{k}^dag = "+" ".join(terms))

# ---- rotations 6.3.36 (spherical) via D^1 ----
print("\n6.3.36 rotated basis (spherical) vs D^1 columns")
import sys,os; sys.path.insert(0,os.path.dirname(__file__))
from wigner_d import wigner_d, beta as B
def D1(a,b,g):
    M=mp.zeros(3,3); order=[1,0,-1]
    for i,mm in enumerate(order):
        for j,mp2 in enumerate(order):
            dd=complex(wigner_d(1,mm,mp2).subs(B,mp.mpf(b)).evalf(25))
            M[i,j]=mp.e**(-1j*mm*a)*mp.mpf(dd.real)*mp.e**(-1j*mp2*g)
    return M
a,b,g=mp.mpf('0.7'),mp.mpf('1.1'),mp.mpf('0.5')
Dm=D1(a,b,g)
cb=mp.cos(b); sb=mp.sin(b)
book_11=mp.matrix([(1+cb)/2*mp.e**(-1j*(a+g)), sb/r2*mp.e**(-1j*g), (1-cb)/2*mp.e**(1j*(a-g))])
book_10=mp.matrix([-sb/r2*mp.e**(-1j*a), cb, sb/r2*mp.e**(1j*a)])
book_1m1=mp.matrix([(1-cb)/2*mp.e**(-1j*(a-g)), -sb/r2*mp.e**(1j*g), (1+cb)/2*mp.e**(1j*(a+g))])
# chi'_{1m'} = sum_m D_{m m'} chi_1m  -> column m' of D
ok&=report("6.3.36 chi'_11", mnorm(mp.matrix([Dm[0,0],Dm[1,0],Dm[2,0]])-book_11))
ok&=report("6.3.36 chi'_10", mnorm(mp.matrix([Dm[0,1],Dm[1,1],Dm[2,1]])-book_10))
ok&=report("6.3.36 chi'_1-1", mnorm(mp.matrix([Dm[0,2],Dm[1,2],Dm[2,2]])-book_1m1))

# ---- helicity 6.3.49 (spherical) : chi_1lam = column lam of D1(phi,theta,0) ----
print("\n6.3.49 helicity (spherical) vs D^1(phi,theta,0)")
th,ph=mp.mpf('1.0'),mp.mpf('0.7')
Dh=D1(ph,th,0); ct=mp.cos(th); st=mp.sin(th)
h11=mp.matrix([(1+ct)/2*mp.e**(-1j*ph), st/r2, (1-ct)/2*mp.e**(1j*ph)])
h10=mp.matrix([-st/r2*mp.e**(-1j*ph), ct, st/r2*mp.e**(1j*ph)])
h1m1=mp.matrix([(1-ct)/2*mp.e**(-1j*ph), -st/r2, (1+ct)/2*mp.e**(1j*ph)])
ok&=report("6.3.49 chi_11(th,ph)", mnorm(mp.matrix([Dh[0,0],Dh[1,0],Dh[2,0]])-h11))
ok&=report("6.3.49 chi_10(th,ph)", mnorm(mp.matrix([Dh[0,1],Dh[1,1],Dh[2,1]])-h10))
ok&=report("6.3.49 chi_1-1(th,ph)", mnorm(mp.matrix([Dh[0,2],Dh[1,2],Dh[2,2]])-h1m1))
# 6.3.50 chi_i(th,ph) from 6.3.44
hx=(mp.matrix([Dh[0,2],Dh[1,2],Dh[2,2]])-mp.matrix([Dh[0,0],Dh[1,0],Dh[2,0]]))/r2
hy=1j*(mp.matrix([Dh[0,2],Dh[1,2],Dh[2,2]])+mp.matrix([Dh[0,0],Dh[1,0],Dh[2,0]]))/r2
hz=mp.matrix([Dh[0,1],Dh[1,1],Dh[2,1]])
book_hx=mp.matrix([-ct/r2*mp.e**(-1j*ph), -st, ct/r2*mp.e**(1j*ph)])
book_hy=mp.matrix([1j/r2*mp.e**(-1j*ph), 0, 1j/r2*mp.e**(1j*ph)])
book_hz=mp.matrix([-st/r2*mp.e**(-1j*ph), ct, st/r2*mp.e**(1j*ph)])
ok&=report("6.3.50 chi_x(th,ph)", mnorm(hx-book_hx))
ok&=report("6.3.50 chi_y(th,ph)", mnorm(hy-book_hy))
ok&=report("6.3.50 chi_z(th,ph)", mnorm(hz-book_hz))

# ---- 6.3.62 <T_2M> in helicity state = (-1)^{1+lam}/(1+|lam|) sqrt(8pi/15) Y_2M ----
print("\n6.3.62 <lam|T_2M|lam> = (-1)^{1+lam}/(1+|lam|) sqrt(8pi/15) Y_2M")
hel={1:mp.matrix([Dh[0,0],Dh[1,0],Dh[2,0]]),0:mp.matrix([Dh[0,1],Dh[1,1],Dh[2,1]]),-1:mp.matrix([Dh[0,2],Dh[1,2],Dh[2,2]])}
w=mp.mpf(0)
for lam in (1,0,-1):
    for M in range(-2,3):
        lhs=(hel[lam].H*Tq[M]*hel[lam])[0]
        rhs=mp.mpf((-1)**(1+lam))/(1+abs(lam))*mp.sqrt(8*mp.pi/15)*mp.spherharm(2,M,th,ph)
        w=max(w,abs(lhs-rhs))
ok&=report("6.3.62 all lam,M", w, tol=mp.mpf('1e-12'))

print("\nRESULT anchors:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
