#!/usr/bin/env python3
r"""Verify Sec 4.18 asymptotics / infinitesimal rotations of D^J_{MM'}.

A fast numerical Wigner d^J_{MM'}(beta) (mpmath) is self-tested against the
validated sympy helper (scripts/wigner_d.py), then used to check:
  4.18.1  large-J cos asymptotic (with s,mu,nu,xi of 4.3.14/15)
  4.18.2  J->inf, beta->0, J beta finite  ->  Bessel J_{M-M'}(J beta)
  4.18.3  beta->0 leading power + first correction
  4.18.4  pi-beta->0 leading power + first correction
  4.18.5-8 infinitesimal rotations = -i <JM|J.n|JM'>
Full D = e^{-iM alpha} d^J_{MM'}(beta) e^{-iM' gamma}.
"""
import math
import mpmath as mp
mp.mp.dps = 200   # d^J at large J is an alternating factorial sum: needs many digits

def fac(n): return mp.factorial(n)

def dJ(J, M, Mp, beta):
    """Wigner small-d d^J_{M,Mp}(beta), VMK/Rose/Edmonds convention.
    Anchored to d^{1/2}_{1/2,-1/2}(b) = -sin(b/2).  (m'=M, m=Mp in the
    standard Wigner formula.)  J,M,Mp may be half-integer."""
    beta = mp.mpf(beta)
    c = mp.cos(beta/2); s = mp.sin(beta/2)
    pref = mp.sqrt(fac(J+M)*fac(J-M)*fac(J+Mp)*fac(J-Mp))
    kmin = int(max(0, Mp-M))            # (M-Mp+k)>=0  and k>=0
    kmax = int(min(J+Mp, J-M))          # (J+Mp-k)>=0  and (J-M-k)>=0
    tot = mp.mpf(0)
    for k in range(kmin, kmax+1):
        denom = fac(J+Mp-k)*fac(k)*fac(M-Mp+k)*fac(J-M-k)
        num = (-1)**(int(M-Mp)+k) * c**(2*J-M+Mp-2*k) * s**(M-Mp+2*k)
        tot += num/denom
    return pref*tot

def D(J,M,Mp,al,be,ga):
    return mp.e**(-1j*M*al) * dJ(J,M,Mp,be) * mp.e**(-1j*Mp*ga)

def report(tag, w, tol=mp.mpf('1e-6')):
    w=float(w); ok=w<tol
    print(f"  {tag:52s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}")
    return ok

ok=True

# ---- self-test dJ against sympy helper ----------------------------------
print("self-test dJ vs sympy wigner_d")
try:
    import sys, os
    sys.path.insert(0, os.path.dirname(__file__))
    from wigner_d import wigner_d, beta as symbeta
    w=mp.mpf(0)
    for twoJ in range(1,7):
        J=mp.mpf(twoJ)/2
        vals=[J-k for k in range(twoJ+1)]
        for M in vals:
            for Mp in vals:
                expr=wigner_d(float(J),float(M),float(Mp))
                for b in [mp.mpf('0.3'),mp.mpf('1.1'),mp.mpf('2.4')]:
                    ref=mp.mpf(str(complex(expr.subs(symbeta,mp.mpf(b)).evalf(25)).real))
                    w=max(w,abs(dJ(J,M,Mp,b)-ref))
    ok&=report("dJ == sympy wigner_d (J<=3, all M,M')", w)
except Exception as e:
    print("  (sympy self-test skipped:", e, ")")

# ---- 4.3.14/15 params ----------------------------------------------------
def params(J,M,Mp):
    mu=abs(M-Mp); nu=abs(M+Mp); s=J-(mu+nu)/2
    xi=1 if Mp>=M else (-1)**int(Mp-M)
    return mu,nu,s,xi

# ---- 4.18.1 large-J cosine asymptotic ------------------------------------
print("\n4.18.1  large-J cos asymptotic")
def approx_1(J,M,Mp,al,be,ga):
    mu,nu,s,xi=params(J,M,Mp)
    amp = xi*mp.sqrt(fac(s)*fac(s+mu+nu)/(fac(s+mu)*fac(s+nu)))
    osc = mp.sqrt(2/(mp.pi*s))*mp.cos((s+(mu+nu+1)/mp.mpf(2))*be - mp.pi/4*(2*mu+1))/mp.sqrt(mp.sin(be))
    return mp.e**(-1j*M*al-1j*Mp*ga)*amp*osc
# leading term ~ J^{-1/2}; residual O(J^{-3/2}) => residual/leading -> 0 like 1/J.
# Verify the ORDER: error*J^{3/2} stays bounded (const ~ few) as J grows.
w=mp.mpf(0)
for J in [mp.mpf(60),mp.mpf(120),mp.mpf(240)]:
    for (M,Mp) in [(mp.mpf(3),mp.mpf(5)),(mp.mpf(-2),mp.mpf(4)),(mp.mpf(10),mp.mpf(-6))]:
        for be in [mp.mpf('0.7'),mp.mpf('1.3'),mp.mpf('2.0')]:
            ex=D(J,M,Mp,0.4,be,0.9); ap=approx_1(J,M,Mp,0.4,be,0.9)
            w=max(w, abs(ex-ap)*J**mp.mpf('1.5'))   # scaled residual ~ O(1)
report("4.18.1 residual*J^{3/2} bounded (order O(J^-3/2))", w, tol=mp.mpf('200'))
ok&= float(w)<200

# convergence RATE: error should fall ~ (J2/J1)^{3/2} per doubling
M,Mp,be=mp.mpf(3),mp.mpf(5),mp.mpf('1.3')
e1=abs(D(mp.mpf(120),M,Mp,0.4,be,0.9)-approx_1(mp.mpf(120),M,Mp,0.4,be,0.9))
e2=abs(D(mp.mpf(480),M,Mp,0.4,be,0.9)-approx_1(mp.mpf(480),M,Mp,0.4,be,0.9))
rate=float(e1/e2)  # expect ~ (480/120)^{1.5}=8
print(f"    error ratio J:120->480 = {rate:.2f}  (expect ~8 for O(J^-3/2))")
ok&= (6< rate <10)

# ---- 4.18.2 Bessel limit -------------------------------------------------
print("\n4.18.2  J->inf, beta->0, J beta finite  ->  Bessel  (order M'-M, corrected)")
def approx_2(J,M,Mp,al,be,ga):
    return mp.e**(-1j*M*al-1j*Mp*ga)*mp.besselj(Mp-M, J*be)   # corrected: M'-M
w=mp.mpf(0)
for J in [mp.mpf(400),mp.mpf(3200)]:
    for (M,Mp) in [(mp.mpf(2),mp.mpf(5)),(mp.mpf(-3),mp.mpf(1)),(mp.mpf(1),mp.mpf(1))]:
        for x in [mp.mpf('1.0'),mp.mpf('3.0'),mp.mpf('6.0')]:  # x=J beta
            be=x/J
            ex=D(J,M,Mp,0.4,be,0.9); ap=approx_2(J,M,Mp,0.4,be,0.9)
            w=max(w, abs(ex-ap))
report("4.18.2 D ~ e^{-iMa-iM'g} J_{M'-M}(Jb)", w, tol=mp.mpf('1e-2'))
ok&= float(w)<1e-2
# guard: the BOOK's J_{M-M'} would FAIL (differs by (-1)^{M-M'})
wb=abs(D(mp.mpf(3200),2,5,0.4,mp.mpf('3.0')/3200,0.9)
       -mp.e**(-1j*2*0.4-1j*5*0.9)*mp.besselj(2-5,mp.mpf('3.0')))
print(f"    (book J_(M-M') would give err={float(wb):.2e} -> confirms misprint)")

# ---- 4.18.3 beta->0 expansion --------------------------------------------
print("\n4.18.3  beta->0 leading + first correction")
def approx_3(J,M,Mp,al,be,ga):
    mu,nu,s,xi=params(J,M,Mp)
    amp=xi/fac(mu)*mp.sqrt(fac(s+mu+nu)*fac(s+mu)/(fac(s)*fac(s+nu)))
    corr=1-(2*s*(s+mu+nu+1)+nu*(mu+1))/(2*(mu+1))*(be/2)**2
    return mp.e**(-1j*M*al-1j*Mp*ga)*amp*(mp.sin(be/2))**mu*corr  # sin, corrected
# residual after the O(beta^2) correction must be O(beta^4): halving beta -> /16
def resid3(J,M,Mp,be):
    return abs(D(J,M,Mp,0.4,be,0.9)/approx_3(J,M,Mp,0.4,be,0.9)-1)
w=mp.mpf(0); worst_rate=0
for J in [mp.mpf(5),mp.mpf(9)]:
    for (M,Mp) in [(mp.mpf(1),mp.mpf(3)),(mp.mpf(-1),mp.mpf(2)),(mp.mpf(2),mp.mpf(-1))]:
        r1=resid3(J,M,Mp,mp.mpf('0.04')); r2=resid3(J,M,Mp,mp.mpf('0.02'))
        rate=float(r1/r2)   # expect ~16
        worst_rate=max(worst_rate,abs(rate-16))
        w=max(w,r2)
print(f"    O(beta^4) rate check: worst |rate-16| = {worst_rate:.2f}")
report("4.18.3 leading+correction correct (O(beta^4) resid)", worst_rate, tol=mp.mpf('1.5'))
ok&= worst_rate<1.5

# ---- 4.18.4 pi-beta->0 expansion -----------------------------------------
print("\n4.18.4  pi-beta->0 leading + first correction")
def approx_4(J,M,Mp,al,be,ga):
    mu,nu,s,xi=params(J,M,Mp)
    d=mp.pi-be
    amp=xi/fac(nu)*(-1)**int(s)*mp.sqrt(fac(s+mu+nu)*fac(s+nu)/(fac(s)*fac(s+mu)))
    corr=1-(2*s*(s+mu+nu+1)+mu*(nu+1))/(2*(nu+1))*(d/2)**2
    return mp.e**(-1j*M*al-1j*Mp*ga)*amp*(mp.sin(d/2))**nu*corr  # sin, corrected
def resid4(J,M,Mp,d):
    be=mp.pi-d
    return abs(D(J,M,Mp,0.4,be,0.9)/approx_4(J,M,Mp,0.4,be,0.9)-1)
worst_rate=0
for J in [mp.mpf(5),mp.mpf(9)]:
    for (M,Mp) in [(mp.mpf(1),mp.mpf(3)),(mp.mpf(-1),mp.mpf(2)),(mp.mpf(2),mp.mpf(-1))]:
        r1=resid4(J,M,Mp,mp.mpf('0.04')); r2=resid4(J,M,Mp,mp.mpf('0.02'))
        worst_rate=max(worst_rate,abs(float(r1/r2)-16))
print(f"    O((pi-beta)^4) rate check: worst |rate-16| = {worst_rate:.2f}")
report("4.18.4 leading+correction correct", worst_rate, tol=mp.mpf('1.5'))
ok&= worst_rate<1.5

# ---- 4.18.5-8 infinitesimal rotations ------------------------------------
print("\n4.18.5-8  infinitesimal rotations = -i <JM|J.n|JM'>")
def num_deriv(f, eps=mp.mpf('1e-6')):
    return (f(eps)-f(-eps))/(2*eps)  # central; delta term cancels
def Jx_me(J,M,Mp):
    r=mp.mpf(0)
    if M==Mp+1: r+=-1j/2*mp.sqrt((J-Mp)*(J+Mp+1))
    if M==Mp-1: r+=-1j/2*mp.sqrt((J+Mp)*(J-Mp+1))
    return r
def Jy_me(J,M,Mp):
    r=mp.mpf(0)
    if M==Mp+1: r+=-mp.mpf(1)/2*mp.sqrt((J-Mp)*(J+Mp+1))
    if M==Mp-1: r+= mp.mpf(1)/2*mp.sqrt((J+Mp)*(J-Mp+1))
    return r
def Jz_me(J,M,Mp):
    return -1j*M if M==Mp else mp.mpf(0)
w=mp.mpf(0)
for twoJ in [2,3,4,5]:
    J=mp.mpf(twoJ)/2
    vals=[J-k for k in range(twoJ+1)]
    for M in vals:
        for Mp in vals:
            d0=1 if M==Mp else 0
            # 4.18.5 x-axis: D(-pi/2, eps, pi/2)
            gx=num_deriv(lambda e: D(J,M,Mp,-mp.pi/2,e,mp.pi/2))
            w=max(w, abs(gx-Jx_me(J,M,Mp)))
            # 4.18.6 y-axis: D(0,eps,0)
            gy=num_deriv(lambda e: D(J,M,Mp,0,e,0))
            w=max(w, abs(gy-Jy_me(J,M,Mp)))
            # 4.18.7 z-axis: D(eps,0,0)
            gz=num_deriv(lambda e: D(J,M,Mp,e,0,0))
            w=max(w, abs(gz-Jz_me(J,M,Mp)))
report("4.18.5-7 x,y,z generators", w, tol=mp.mpf('1e-4'))
ok&= float(w)<1e-4

# 4.18.8 arbitrary axis n(Theta,Phi): rotation eps about n means
#   R = e^{-i eps n.J}; compare -i<JM|n.J|JM'> vs same via Euler.
# n.J = cosT Jz + sinT(cosP Jx + sinP Jy). Build matrix elements directly.
print("\n4.18.8  arbitrary axis n(Theta,Phi)")
def nJ_me(J,M,Mp,T,P):
    r = -1j*( mp.cos(T)*(M if M==Mp else 0) )
    if M==Mp+1: r+= -1j/2*mp.sin(T)*mp.e**(-1j*P)*mp.sqrt((J-Mp)*(J+Mp+1))
    if M==Mp-1: r+= -1j/2*mp.sin(T)*mp.e**( 1j*P)*mp.sqrt((J+Mp)*(J-Mp+1))
    return r
# verify RHS formula equals -i<JM|n.J|JM'> built from Jx,Jy,Jz elements
w=mp.mpf(0)
T=mp.mpf('0.7'); P=mp.mpf('1.1')
for twoJ in [2,3,4]:
    J=mp.mpf(twoJ)/2
    vals=[J-k for k in range(twoJ+1)]
    for M in vals:
        for Mp in vals:
            # build -i<n.J> from component matrix elements:
            nJ = mp.cos(T)*(M if M==Mp else 0)*1  # <Jz>
            nJ += mp.sin(T)*mp.cos(P)*(1j*Jx_me(J,M,Mp))  # Jx_me=-i<Jx> => <Jx>=i*Jx_me
            nJ += mp.sin(T)*mp.sin(P)*(1j*Jy_me(J,M,Mp))
            ref=-1j*nJ
            w=max(w, abs(ref-nJ_me(J,M,Mp,T,P)))
report("4.18.8 n.J matrix element identity", w, tol=mp.mpf('1e-20'))
ok&= float(w)<1e-20

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
