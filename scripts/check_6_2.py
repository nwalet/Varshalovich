#!/usr/bin/env python3
r"""Verify Sec 6.2 (spin functions for S=1/2) as explicit 2x2 linear algebra.

Conventions (from the section itself):
  basis alpha=chi_{1/2,1/2}=(1,0)^T, beta=chi_{1/2,-1/2}=(0,1)^T
  S_i = sigma_i/2 ; spherical S_0=S_z, S_{+1}=-(S_x+iS_y)/sqrt2, S_{-1}=(S_x-iS_y)/sqrt2
  D^{1/2}_{mm'}(a,b,g) = e^{-i m a} d^{1/2}_{mm'}(b) e^{-i m' g},
     d^{1/2}=[[cos b/2, -sin b/2],[sin b/2, cos b/2]]  (rows/cols m=+1/2,-1/2)
  U^{1/2}(w;T,P) = exp(-i w n.S),  n=(sinT cosP, sinT sinP, cosT)
  clebsch{a}{b}{c}{d}{e}{f} = C^{e f}_{a b, c d} = clebsch_gordan(a,c,e,b,d,f)
Checks the book's explicit forms/relations 6.2.9-6.2.31 against these.
"""
import mpmath as mp
mp.mp.dps = 30
I2 = mp.matrix([[1,0],[0,1]])
sx = mp.matrix([[0,1],[1,0]])/2
sy = mp.matrix([[0,-1j],[1j,0]])/2
sz = mp.matrix([[1,0],[0,-1]])/2
Sph = {0: sz, 1: -(sx+1j*sy)/mp.sqrt(2), -1: (sx-1j*sy)/mp.sqrt(2)}
ket = {mp.mpf(1)/2: mp.matrix([1,0]), mp.mpf(-1)/2: mp.matrix([0,1])}
def dag(v):  # conjugate transpose of column or matrix
    return v.H
def cg(a,b,c,d,e,f):
    from sympy.physics.quantum.cg import CG
    from sympy import Rational as R, N
    from fractions import Fraction
    def rr(x): return R(Fraction(float(x)).limit_denominator(8))
    return complex(N(CG(rr(a),rr(b),rr(c),rr(d),rr(e),rr(f)).doit()))

def d12(b):
    c=mp.cos(b/2); s=mp.sin(b/2)
    return mp.matrix([[c,-s],[s,c]])
def D12(a,b,g):
    d=d12(b); ph=[mp.e**(-1j*mp.mpf(1)/2*a), mp.e**(1j*mp.mpf(1)/2*a)]
    phg=[mp.e**(-1j*mp.mpf(1)/2*g), mp.e**(1j*mp.mpf(1)/2*g)]
    return mp.matrix([[ph[i]*d[i,j]*phg[j] for j in range(2)] for i in range(2)])
def U12(w,T,P):
    nx=mp.sin(T)*mp.cos(P); ny=mp.sin(T)*mp.sin(P); nz=mp.cos(T)
    return mp.expm(-1j*w*(nx*sx+ny*sy+nz*sz))

def mnorm(M):
    return max(abs(M[i,j]) for i in range(M.rows) for j in range(M.cols))
def report(tag, w, tol=mp.mpf('1e-24')):
    w=float(w); ok=w<tol
    print(f"  {tag:46s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}")
    return ok
ok=True
half=mp.mpf(1)/2

# ---- 6.2.9/10 product expansions ----
print("6.2.10 product expansions chi chi^dagger")
aa=ket[half]*dag(ket[half]); ab=ket[half]*dag(ket[-half])
ba=ket[-half]*dag(ket[half]); bb=ket[-half]*dag(ket[-half])
ok&=report("aa^d = I/2 + S_z", mnorm(aa-(I2/2+sz)))
ok&=report("ab^d = S_x+iS_y = -sqrt2 S_+1", mnorm(ab-(sx+1j*sy))+mnorm(ab+mp.sqrt(2)*Sph[1]))
ok&=report("ba^d = S_x-iS_y = sqrt2 S_-1", mnorm(ba-(sx-1j*sy))+mnorm(ba-mp.sqrt(2)*Sph[-1]))
ok&=report("bb^d = I/2 - S_z", mnorm(bb-(I2/2-sz)))
# 6.2.9 general: chi_m chi_m'^d = 1/2 delta I - sqrt3 C^{1/2 m}_{1 mu,1/2 m'} S_mu
w=mp.mpf(0)
for m in [half,-half]:
    for mp_ in [half,-half]:
        lhs=ket[m]*dag(ket[mp_])
        rhs=half*(1 if m==mp_ else 0)*I2
        for mu in [1,0,-1]:
            rhs=rhs - mp.sqrt(3)*cg(1,mu,half,mp_,half,m)*Sph[mu]
        w=max(w,mnorm(lhs-rhs))
ok&=report("6.2.9 general CG form", w, tol=mp.mpf('1e-14'))

# ---- 6.2.11 cartesian action ----
print("\n6.2.11/13 spin-operator action on basis")
w=mp.mpf(0)
w=max(w, mnorm(sx*ket[half]-half*ket[-half]), mnorm(sx*ket[-half]-half*ket[half]))
w=max(w, mnorm(sy*ket[half]-1j/2*ket[-half]), mnorm(sy*ket[-half]+1j/2*ket[half]))
w=max(w, mnorm(sz*ket[half]-half*ket[half]), mnorm(sz*ket[-half]+half*ket[-half]))
ok&=report("6.2.11 S_x,S_y,S_z chi", w)
# 6.2.13 spherical action
w=mp.mpf(0)
w=max(w, mnorm(Sph[1]*ket[half]), mnorm(Sph[1]*ket[-half]+1/mp.sqrt(2)*ket[half]))
w=max(w, mnorm(Sph[0]*ket[half]-half*ket[half]), mnorm(Sph[0]*ket[-half]+half*ket[-half]))
w=max(w, mnorm(Sph[-1]*ket[half]-1/mp.sqrt(2)*ket[-half]), mnorm(Sph[-1]*ket[-half]))
ok&=report("6.2.13 S_+1,S_0,S_-1 chi", w)
# 6.2.12 general: S_mu chi_m = -sqrt3/2 C^{1/2 m'}_{1 mu,1/2 m} chi_m'
w=mp.mpf(0)
for mu in [1,0,-1]:
    for m in [half,-half]:
        lhs=Sph[mu]*ket[m]
        rhs=mp.matrix([0,0])
        for mpp in [half,-half]:
            rhs=rhs - mp.sqrt(3)/2*cg(1,mu,half,m,half,mpp)*ket[mpp]
        w=max(w,mnorm(lhs-rhs))
ok&=report("6.2.12 general CG form", w, tol=mp.mpf('1e-14'))

# ---- 6.2.17 rotated Euler basis (gives correct exponents) ----
print("\n6.2.17 rotated basis (Euler)  vs  D^{1/2} columns")
a,b,g=mp.mpf('0.7'),mp.mpf('1.1'),mp.mpf('0.5')
D=D12(a,b,g)
# chi'_{1/2 m'} = sum_m chi_m D_{m m'}: column m' of D
chip_half=mp.matrix([D[0,0],D[1,0]])   # m'=+1/2
chip_mhalf=mp.matrix([D[0,1],D[1,1]])  # m'=-1/2
book_half=mp.matrix([mp.cos(b/2)*mp.e**(-1j*(a+g)/2), mp.sin(b/2)*mp.e**(1j*(a-g)/2)])
book_mhalf=mp.matrix([-mp.sin(b/2)*mp.e**(-1j*(a-g)/2), mp.cos(b/2)*mp.e**(1j*(a+g)/2)])
ok&=report("6.2.17 chi'_+ (corrected exponents)", mnorm(chip_half-book_half))
ok&=report("6.2.17 chi'_- (corrected exponents)", mnorm(chip_mhalf-book_mhalf))

# ---- 6.2.18 rotated basis (axis-angle) ----
print("\n6.2.18 rotated basis (omega;Theta,Phi)  vs  U^{1/2} columns")
w_,T,P=mp.mpf('0.9'),mp.mpf('1.3'),mp.mpf('0.6')
U=U12(w_,T,P)
chip_h=mp.matrix([U[0,0],U[1,0]]); chip_mh=mp.matrix([U[0,1],U[1,1]])
book_h=mp.matrix([mp.cos(w_/2)-1j*mp.sin(w_/2)*mp.cos(T), -1j*mp.sin(w_/2)*mp.sin(T)*mp.e**(1j*P)])
book_mh=mp.matrix([-1j*mp.sin(w_/2)*mp.sin(T)*mp.e**(-1j*P), mp.cos(w_/2)+1j*mp.sin(w_/2)*mp.cos(T)])
ok&=report("6.2.18 chi'_+ ", mnorm(chip_h-book_h))
ok&=report("6.2.18 chi'_- ", mnorm(chip_mh-book_mh))

# ---- 6.2.22 helicity functions ----
print("\n6.2.22 helicity functions  vs  sum_m D^{1/2}_{m lam}(phi,theta,0) chi_m")
th,ph=mp.mpf('1.0'),mp.mpf('0.7')
Dh=D12(ph,th,0)
hel_p=mp.matrix([Dh[0,0],Dh[1,0]])   # lambda=+1/2
hel_m=mp.matrix([Dh[0,1],Dh[1,1]])   # lambda=-1/2
book_hp=mp.matrix([mp.cos(th/2)*mp.e**(-1j*ph/2), mp.sin(th/2)*mp.e**(1j*ph/2)])
book_hm=mp.matrix([-mp.sin(th/2)*mp.e**(-1j*ph/2), mp.cos(th/2)*mp.e**(1j*ph/2)])
ok&=report("6.2.22 chi_+1/2(theta,phi)", mnorm(hel_p-book_hp))
ok&=report("6.2.22 chi_-1/2(theta,phi)", mnorm(hel_m-book_hm))

# ---- 6.2.24/25 helicity orthonormality/completeness ----
print("\n6.2.24/25 helicity orthonormality & completeness")
ok&=report("6.2.24 <lam'|lam>=delta", abs((dag(hel_p)*hel_m)[0])+abs((dag(hel_p)*hel_p)[0]-1))
comp=hel_p*dag(hel_p)+hel_m*dag(hel_m)
ok&=report("6.2.25 sum |lam><lam| = I", mnorm(comp-I2))

# ---- 6.2.27 helicity product matrices ----
print("\n6.2.27 helicity products (explicit 2x2)")
nx=mp.sin(th)*mp.cos(ph); ny=mp.sin(th)*mp.sin(ph); nz=mp.cos(th)
nS=nx*sx+ny*sy+nz*sz
pp=hel_p*dag(hel_p); pm=hel_p*dag(hel_m); mp_=hel_m*dag(hel_p); mm=hel_m*dag(hel_m)
ok&=report("6.2.27 |+><+| = I/2 + n.S", mnorm(pp-(I2/2+nS)))
ok&=report("6.2.27 |-><-| = I/2 - n.S", mnorm(mm-(I2/2-nS)))
book_pp=mp.matrix([[1+mp.cos(th), mp.sin(th)*mp.e**(-1j*ph)],[mp.sin(th)*mp.e**(1j*ph),1-mp.cos(th)]])/2
ok&=report("6.2.27 |+><+| explicit matrix", mnorm(pp-book_pp))
book_pm=mp.matrix([[-mp.sin(th),(mp.cos(th)+1)*mp.e**(-1j*ph)],[(mp.cos(th)-1)*mp.e**(1j*ph),mp.sin(th)]])/2
ok&=report("6.2.27 |+><-| explicit matrix", mnorm(pm-book_pm))
book_mp=mp.matrix([[-mp.sin(th),(mp.cos(th)-1)*mp.e**(-1j*ph)],[(mp.cos(th)+1)*mp.e**(1j*ph),mp.sin(th)]])/2
ok&=report("6.2.27 |-><+| explicit matrix", mnorm(mp_-book_mp))

# ---- 6.2.29 helicity matrix elements of spherical S_mu ----
print("\n6.2.29 <lam'|S_mu|lam> (helicity)")
hel={half:hel_p,-half:hel_m}
def me(lp,mu,lam): return (dag(hel[lp])*Sph[mu]*hel[lam])[0]
book29={
 (half,1,half): -mp.sin(th)/(2*mp.sqrt(2))*mp.e**(1j*ph),
 (half,1,-half): -(1+mp.cos(th))/(2*mp.sqrt(2))*mp.e**(1j*ph),
 (half,0,half): half*mp.cos(th),
 (half,0,-half): -half*mp.sin(th),
 (half,-1,half): mp.sin(th)/(2*mp.sqrt(2))*mp.e**(-1j*ph),
 (half,-1,-half): -(1-mp.cos(th))/(2*mp.sqrt(2))*mp.e**(-1j*ph),
 (-half,1,half): (1-mp.cos(th))/(2*mp.sqrt(2))*mp.e**(1j*ph),
 (-half,1,-half): mp.sin(th)/(2*mp.sqrt(2))*mp.e**(1j*ph),
 (-half,0,half): -half*mp.sin(th),
 (-half,0,-half): -half*mp.cos(th),
 (-half,-1,half): (1+mp.cos(th))/(2*mp.sqrt(2))*mp.e**(-1j*ph),
 (-half,-1,-half): -mp.sin(th)/(2*mp.sqrt(2))*mp.e**(-1j*ph),
}
w=mp.mpf(0)
for k,v in book29.items():
    w=max(w, abs(me(k[0],k[1],k[2])-v))
ok&=report("6.2.29 all 12 helicity S_mu elements", w)

# ---- 6.2.30 cartesian helicity matrix elements ----
print("\n6.2.30 <lam'|S_i|lam> (helicity, cartesian)")
def meC(lp,Si,lam): return (dag(hel[lp])*Si*hel[lam])[0]
Smap={'x':sx,'y':sy,'z':sz}
book30=[
 ('+','x','+', half*mp.sin(th)*mp.cos(ph)),
 ('+','x','-', half*(mp.cos(th)*mp.cos(ph)+1j*mp.sin(ph))),
 ('+','y','+', half*mp.sin(th)*mp.sin(ph)),
 ('+','y','-', half*(mp.cos(th)*mp.sin(ph)-1j*mp.cos(ph))),
 ('+','z','+', half*mp.cos(th)),
 ('+','z','-', -half*mp.sin(th)),
 ('-','x','+', half*(mp.cos(th)*mp.cos(ph)-1j*mp.sin(ph))),
 ('-','x','-', -half*mp.sin(th)*mp.cos(ph)),
 ('-','y','+', half*(mp.cos(th)*mp.sin(ph)+1j*mp.cos(ph))),
 ('-','y','-', -half*mp.sin(th)*mp.sin(ph)),
 ('-','z','+', -half*mp.sin(th)),
 ('-','z','-', -half*mp.cos(th)),
]
lam={'+':half,'-':-half}
w=mp.mpf(0)
for lp,Si,lm,v in book30:
    w=max(w, abs(meC(lam[lp],Smap[Si],lam[lm])-v))
ok&=report("6.2.30 all 12 helicity cartesian elements", w)

# ---- 6.2.31 diagonal <S> = lambda n ----
print("\n6.2.31 diagonal <lam|S|lam> = lambda n")
w=mp.mpf(0)
nvec=mp.matrix([nx,ny,nz])
for lam_,sgn in [(half,1),(-half,-1)]:
    Svec=mp.matrix([(dag(hel[lam_])*Si*hel[lam_])[0] for Si in (sx,sy,sz)])
    w=max(w, mnorm(Svec-lam_*nvec))
ok&=report("6.2.31 <S>=lambda n", w)

# ---- 6.2.39/41 spin direction n from a generic spinor ----
print("\n6.2.39/41 n-vector from spinor a=(a^1/2,a^-1/2)")
aU=mp.mpf('0.6')+0.3j; bU=mp.mpf('0.5')-0.4j
nrm=mp.sqrt(abs(aU)**2+abs(bU)**2); aU/=nrm; bU/=nrm
chi=mp.matrix([aU,bU])
# n = <sigma> = chi^dag sigma chi  (P=2<S>=<sigma>)
nx=(dag(chi)*(2*sx)*chi)[0]; ny=(dag(chi)*(2*sy)*chi)[0]; nz=(dag(chi)*(2*sz)*chi)[0]
ok&=report("6.2.39 n_x=2Re(a*b)", abs(nx-2*(mp.conj(aU)*bU).real))
ok&=report("6.2.39 n_y=2Im(a*b)", abs(ny-2*(mp.conj(aU)*bU).imag))
ok&=report("6.2.39 n_z=|a|^2-|b|^2", abs(nz-(abs(aU)**2-abs(bU)**2)))
# spherical 6.2.41
np1=-(nx+1j*ny)/mp.sqrt(2); n0=nz; nm1=(nx-1j*ny)/mp.sqrt(2)
ok&=report("6.2.41 n_+1=-sqrt2 a* b", abs(np1-(-mp.sqrt(2)*mp.conj(aU)*bU)))
ok&=report("6.2.41 n_-1=sqrt2 b* a", abs(nm1-(mp.sqrt(2)*mp.conj(bU)*aU)))
# 6.2.40 general CG form: n_mu = sqrt3 sum C^{1/2 m'}_{1/2 m,1 mu} a^{m'*} a^{m}
kets={half:aU,-half:bU}
for mu,ref in [(1,np1),(0,n0),(-1,nm1)]:
    val=mp.mpf(0)
    for m in [half,-half]:
        for mpp in [half,-half]:
            val+=mp.sqrt(3)*cg(half,m,1,mu,half,mpp)*mp.conj(kets[mpp])*kets[m]
    ok&=report(f"6.2.40 n_{mu} CG form", abs(val-ref), tol=mp.mpf('1e-14'))

# ---- 6.2.50-52 polarization vector from density matrix ----
print("\n6.2.50-52 polarization from rho")
rho=chi*dag(chi)   # pure state
def rr(s,sp):  # rho_{s,s'} with index +1/2->0, -1/2->1
    idx={half:0,-half:1}; return rho[idx[s],idx[sp]]
Px=(rr(half,-half)+rr(-half,half)); Py=1j*(rr(half,-half)-rr(-half,half)); Pz=rr(half,half)-rr(-half,-half)
ok&=report("6.2.52 P_x", abs(Px-nx)); ok&=report("6.2.52 P_y", abs(Py-ny)); ok&=report("6.2.52 P_z", abs(Pz-nz))
Pp1=-mp.sqrt(2)*rr(-half,half); P0=rr(half,half)-rr(-half,-half); Pm1=mp.sqrt(2)*rr(half,-half)
ok&=report("6.2.51 P_+1=-sqrt2 rho_-1/2,1/2", abs(Pp1-np1))
ok&=report("6.2.51 P_-1=sqrt2 rho_1/2,-1/2", abs(Pm1-nm1))
# 6.2.50 CG form
for mu,ref in [(1,Pp1),(0,P0),(-1,Pm1)]:
    val=mp.mpf(0)
    for s in [half,-half]:
        for sp in [half,-half]:
            val+=mp.sqrt(3)*cg(half,s,1,mu,half,sp)*rr(s,sp)
    ok&=report(f"6.2.50 P_{mu} CG form", abs(val-ref), tol=mp.mpf('1e-14'))

# ---- 6.2.46 contravariant transform (does the inverse need a conjugate?) ----
print("\n6.2.46 contravariant spinor transform")
a2,b2,g2=mp.mpf('0.7'),mp.mpf('1.1'),mp.mpf('0.5')
Dm=D12(a2,b2,g2)
avec=mp.matrix([aU,bU])   # a^n, n=+1/2,-1/2
# forward:  a'^m = sum_n D_{mn} a^n
aprime=Dm*avec
# book's 2nd line as printed:  a^n = sum_m D_{mn} a'^m  (NO conjugate)
back_noconj=mp.matrix([sum(Dm[m,n]*aprime[m] for m in range(2)) for n in range(2)])
# correct inverse (unitary): a^n = sum_m D*_{mn} a'^m
back_conj=mp.matrix([sum(mp.conj(Dm[m,n])*aprime[m] for m in range(2)) for n in range(2)])
print(f"    printed (no conj) recovers a? worst={float(mnorm(back_noconj-avec)):.2e}")
print(f"    with conjugate     recovers a? worst={float(mnorm(back_conj-avec)):.2e}")
ok&=report("6.2.46 inverse needs D* (conjugate)", mnorm(back_conj-avec))

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
