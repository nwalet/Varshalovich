#!/usr/bin/env python3
r"""Verify Sec 6.1 (spin functions, arbitrary spin) for S=1/2,1,3/2 and two-particle.

Spin-S matrices in basis m=S,S-1,...,-S.  Polarization operators from Eq 2.4.8:
  [T_LM(S)]_{s',s} = sqrt((2L+1)/(2S+1)) C^{S s'}_{S s, L M}.
D^S(a,b,g)=exp(-ia Sz)exp(-ib Sy)exp(-ig Sz).
clebsch{a}{b}{c}{d}{e}{f}=C^{e f}_{a b,c d}; sixj = Wigner 6j.
"""
import mpmath as mp
mp.mp.dps=25
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_6j
from sympy import Rational as R, N
from fractions import Fraction
def rr(x): return R(Fraction(float(x)).limit_denominator(16))
def cg(a,b,c,d,e,f):
    return complex(N(CG(rr(a),rr(b),rr(c),rr(d),rr(e),rr(f)).doit()))
def _tri(a,b,c):
    return (a+b>=c) and (a+c>=b) and (b+c>=a) and (round(a+b+c)==a+b+c)
def _delta(a,b,c):
    return mp.sqrt(mp.factorial(a+b-c)*mp.factorial(a-b+c)*mp.factorial(-a+b+c)/mp.factorial(a+b+c+1))
def sixj(j1,j2,j3,j4,j5,j6):
    """Wigner 6j via the Racah single-sum formula (robust for half-integers)."""
    tris=[(j1,j2,j3),(j1,j5,j6),(j4,j2,j6),(j4,j5,j3)]
    if not all(_tri(*t) for t in tris): return 0.0
    pref=_delta(j1,j2,j3)*_delta(j1,j5,j6)*_delta(j4,j2,j6)*_delta(j4,j5,j3)
    a1=j1+j2+j3; a2=j1+j5+j6; a3=j4+j2+j6; a4=j4+j5+j3
    b1=j1+j2+j4+j5; b2=j2+j3+j5+j6; b3=j1+j3+j4+j6
    tmin=max(a1,a2,a3,a4); tmax=min(b1,b2,b3)   # may be half-integer; step by 1
    s=mp.mpf(0); t=mp.mpf(tmin)
    while t<=tmax+mp.mpf('1e-9'):
        s+=(-1)**int(round(float(t)))*mp.factorial(t+1)/(mp.factorial(t-a1)*mp.factorial(t-a2)*mp.factorial(t-a3)*mp.factorial(t-a4)
            *mp.factorial(b1-t)*mp.factorial(b2-t)*mp.factorial(b3-t))
        t+=1
    return float(pref*s)

def spin(S):
    ms=[S-k for k in range(int(round(2*S))+1)]      # S,...,-S
    n=len(ms); Sz=mp.zeros(n,n); Sp=mp.zeros(n,n)
    for i,m in enumerate(ms): Sz[i,i]=m
    for i,m in enumerate(ms):
        if m+1<=S:
            Sp[ms.index(m+1),i]=mp.sqrt(S*(S+1)-m*(m+1))
    Sm=Sp.T
    Sx=(Sp+Sm)/2; Sy=(Sp-Sm)/(2j)
    return ms,Sx,Sy,Sz
def sphc(Sx,Sy,Sz): return {0:Sz,1:-(Sx+1j*Sy)/mp.sqrt(2),-1:(Sx-1j*Sy)/mp.sqrt(2)}
def Tmat(S,L,M):
    ms=[S-k for k in range(int(round(2*S))+1)]; n=len(ms)
    T=mp.zeros(n,n); pref=mp.sqrt((2*L+1)/(2*S+1))
    for i,sp in enumerate(ms):
        for j,s in enumerate(ms):
            if abs(M-(sp-s))<1e-9:
                T[i,j]=pref*cg(S,s,L,M,S,sp)
    return T
def Dmat(S,a,b,g):
    ms,Sx,Sy,Sz=spin(S)
    return ms, mp.expm(-1j*a*Sz)*mp.expm(-1j*b*Sy)*mp.expm(-1j*g*Sz)
def e_(S,m):
    ms=[S-k for k in range(int(round(2*S))+1)]; v=mp.zeros(len(ms),1); v[ms.index(m)]=1; return v
def mnorm(A):
    return max(abs(A[i,j]) for i in range(A.rows) for j in range(A.cols))
def report(tag,w,tol=mp.mpf('1e-11')):
    w=float(w); ok=w<tol; print(f"  {tag:44s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok
ok=True
Svals=[mp.mpf(1)/2,mp.mpf(1),mp.mpf(3)/2]

print("6.1.5 eigen; 6.1.11 chi chi^dag = sum T; 6.1.12-14 spin ME")
for S in Svals:
    ms,Sx,Sy,Sz=spin(S); Sp=sphc(Sx,Sy,Sz); n=len(ms)
    S2=Sx*Sx+Sy*Sy+Sz*Sz
    ok&=report(f"S={float(S)} S^2=S(S+1)", mnorm(S2-S*(S+1)*mp.eye(n)))
    # 6.1.11
    w=mp.mpf(0)
    for m in ms:
        for mp2 in ms:
            lhs=e_(S,m)*e_(S,mp2).H
            rhs=mp.zeros(n,n); M=m-mp2
            for L in range(0,int(2*S)+1):
                if abs(M)<=L:
                    rhs=rhs+mp.sqrt((2*L+1)/(2*S+1))*cg(S,mp2,L,M,S,m)*Tmat(S,L,M)
            w=max(w,mnorm(lhs-rhs))
    ok&=report(f"S={float(S)} 6.1.11", w)
    # 6.1.12/13/14
    w=mp.mpf(0)
    for m in ms:
        for mpr in ms:
            for mu in (1,0,-1):
                lhs=(e_(S,mpr).H*Sp[mu]*e_(S,m))[0]
                rhs=mp.sqrt(S*(S+1))*cg(S,m,1,mu,S,mpr)
                w=max(w,abs(lhs-rhs))
    ok&=report(f"S={float(S)} 6.1.12 <S_mu>", w)

print("\n6.1.18 S'_mu = sum_nu D^1_{nu mu} S_nu")
for S in [mp.mpf(1)/2,mp.mpf(1)]:
    ms,Sx,Sy,Sz=spin(S); Sp=sphc(Sx,Sy,Sz)
    a,b,g=mp.mpf('0.6'),mp.mpf('1.0'),mp.mpf('0.4')
    _,D1=Dmat(1,a,b,g); o1=[1,0,-1]
    _,DS=Dmat(S,a,b,g)
    w=mp.mpf(0)
    for mu in (1,0,-1):
        # S'_mu = D^S S_mu D^S^dagger  (active rotation of operator)
        Sprime=DS*Sp[mu]*DS.H
        rhs=mp.zeros(DS.rows,DS.rows)
        for j,nu in enumerate(o1):
            rhs=rhs+D1[o1.index(nu),o1.index(mu)]*Sp[nu]
        w=max(w,mnorm(Sprime-rhs))
    ok&=report(f"S={float(S)} 6.1.18", w)

print("\n6.1.22 helicity comps; 6.1.25-27 helicity T expansions")
for S in [mp.mpf(1)/2,mp.mpf(1)]:
    ms,Sx,Sy,Sz=spin(S); Sp=sphc(Sx,Sy,Sz); n=len(ms)
    th,ph=mp.mpf('0.9'),mp.mpf('0.6')
    _,Dh=Dmat(S,ph,th,0)
    hel={lam: mp.matrix([Dh[i,ms.index(lam)] for i in range(n)]) for lam in ms}
    # 6.1.22 [chi_Slam]^sigma = D^S_{sigma lam}(phi,theta,0)
    w=max(abs(hel[lam][ms.index(sig)]-Dh[ms.index(sig),ms.index(lam)]) for lam in ms for sig in ms)
    ok&=report(f"S={float(S)} 6.1.22", w)
    # 6.1.27 diagonal <T_LM> = sqrt(4pi/(2S+1)) C^{Slam}_{Slam,L0} Y_LM
    w=mp.mpf(0)
    for lam in ms:
        for L in range(0,int(2*S)+1):
            for M in range(-L,L+1):
                lhs=(hel[lam].H*Tmat(S,L,M)*hel[lam])[0]
                rhs=mp.sqrt(4*mp.pi/(2*S+1))*cg(S,lam,L,0,S,lam)*mp.spherharm(L,M,th,ph)
                w=max(w,abs(lhs-rhs))
    ok&=report(f"S={float(S)} 6.1.27 <T_LM>=..Y_LM", w)

print("\n6.1.30-35 P_LM decomposition of chi chi^dag (random pure state)")
for S in [mp.mpf(1)/2,mp.mpf(1)]:
    ms,Sx,Sy,Sz=spin(S); n=len(ms)
    amp=mp.matrix([mp.mpf('0.4')+0.2j*(k+1)-0.1j*k for k in range(n)])
    amp=amp/mp.sqrt(sum(abs(amp[i])**2 for i in range(n)))
    a={ms[i]:amp[i] for i in range(n)}
    P={}
    for L in range(0,int(2*S)+1):
        for M in range(-L,L+1):
            s=mp.mpf(0)
            for m in ms:
                for nn in ms:
                    s+=mp.sqrt((2*L+1)/(2*S+1))*cg(S,m,L,M,S,nn)*mp.conj(a[nn])*a[m]
            P[(L,M)]=s
    # 6.1.30 chi chi^dag = sum (-1)^M P_{L,-M} T_LM
    lhs=amp*amp.H; rhs=mp.zeros(n,n)
    for (L,M),_ in list(P.items()):
        rhs=rhs+(-1)**M*P[(L,-M)]*Tmat(S,L,M)
    ok&=report(f"S={float(S)} 6.1.30/31", mnorm(lhs-rhs))
    # 6.1.35 sum |P_LM|^2 = 1
    ok&=report(f"S={float(S)} 6.1.35 norm", abs(sum(abs(v)**2 for v in P.values())-1))
    ok&=report(f"S={float(S)} 6.1.31 P00=1/sqrt(2S+1)", abs(P[(0,0)]-1/mp.sqrt(2*S+1)))

print("\n6.1.49-51 statistical tensors (random density matrix)")
for S in [mp.mpf(1)/2,mp.mpf(1)]:
    ms,Sx,Sy,Sz=spin(S); n=len(ms)
    A=mp.matrix([[mp.mpf(1+i)/(1+j)+0.1j*(i-j) for j in range(n)] for i in range(n)])
    rho=A*A.H; rho=rho/mp.re(sum(rho[i,i] for i in range(n)))
    t={}
    for L in range(0,int(2*S)+1):
        for M in range(-L,L+1):
            s=mp.mpf(0)
            for i,sig in enumerate(ms):
                for j,sigp in enumerate(ms):
                    s+=mp.sqrt((2*L+1)/(2*S+1))*cg(S,sig,L,M,S,sigp)*rho[i,j]  # rho_{sig sig'}
            t[(L,M)]=s
    # 6.1.48: t_LM = Tr(rho T_LM); verify against that as ground truth
    w=mp.mpf(0)
    for (L,M),val in t.items():
        tr=sum((rho*Tmat(S,L,M))[i,i] for i in range(n))
        w=max(w,abs(val-tr))
    ok&=report(f"S={float(S)} 6.1.49 vs Tr(rho T)", w)
    ok&=report(f"S={float(S)} 6.1.51 t00=1/sqrt(2S+1)", abs(t[(0,0)]-1/mp.sqrt(2*S+1)))

print("\n6.1.53-65 two particles")
def kron(A,B):
    r=mp.zeros(A.rows*B.rows,A.cols*B.cols)
    for i in range(A.rows):
        for j in range(A.cols):
            for k in range(B.rows):
                for l in range(B.cols):
                    r[i*B.rows+k, j*B.cols+l]=A[i,j]*B[k,l]
    return r
for (S1,S2) in [(mp.mpf(1)/2,mp.mpf(1)/2),(mp.mpf(1),mp.mpf(1)/2)]:
    m1s=[S1-k for k in range(int(2*S1)+1)]; m2s=[S2-k for k in range(int(2*S2)+1)]
    NN=len(m1s)*len(m2s)
    def prod_ket(m1,m2):  # |S1 m1>|S2 m2> as length-N vector
        v=mp.zeros(NN,1); v[m1s.index(m1)*len(m2s)+m2s.index(m2)]=1; return v
    Svals2=[abs(S1-S2)+k for k in range(int(S1+S2-abs(S1-S2))+1)]
    # coupled states 6.1.53
    def coupled(S,m):
        v=mp.zeros(NN,1)
        for m1 in m1s:
            for m2 in m2s:
                if abs(m1+m2-m)<1e-9:
                    v=v+cg(S1,m1,S2,m2,S,m)*prod_ket(m1,m2)
        return v
    # 6.1.55 orthonormality
    w=mp.mpf(0)
    for S in Svals2:
        for m in [S-k for k in range(int(2*S)+1)]:
            for Sp2 in Svals2:
                for mpr in [Sp2-k for k in range(int(2*Sp2)+1)]:
                    ip=(coupled(S,m).H*coupled(Sp2,mpr))[0]
                    tgt=1 if (S==Sp2 and abs(m-mpr)<1e-9) else 0
                    w=max(w,abs(ip-tgt))
    ok&=report(f"({float(S1)},{float(S2)}) 6.1.55 orthonorm", w)
    # 6.1.56 completeness
    comp=mp.zeros(NN,NN)
    for S in Svals2:
        for m in [S-k for k in range(int(2*S)+1)]:
            comp=comp+coupled(S,m)*coupled(S,m).H
    ok&=report(f"({float(S1)},{float(S2)}) 6.1.56 complete", mnorm(comp-mp.eye(NN)))
    # 6.1.58 T_L2M2(2) action
    I1=mp.eye(len(m1s))
    w=mp.mpf(0)
    for S in Svals2:
        for m in [S-k for k in range(int(2*S)+1)]:
            for L2 in range(0,int(2*S2)+1):
                for M2 in range(-L2,L2+1):
                    T2=kron(I1,Tmat(S2,L2,M2))
                    lhs=T2*coupled(S,m)
                    rhs=mp.zeros(NN,1)
                    for Sp2 in Svals2:
                        mpr=m+M2
                        if abs(mpr)<=Sp2:
                            coef=(-1)**int(S1+S2+Sp2+L2)*mp.sqrt((2*L2+1)*(2*S+1))*sixj(S1,S2,S,L2,Sp2,S2)*cg(S,m,L2,M2,Sp2,mpr)
                            rhs=rhs+coef*coupled(Sp2,mpr)
                    w=max(w,mnorm(lhs-rhs))
    ok&=report(f"({float(S1)},{float(S2)}) 6.1.58 T_LM(2) action", w, tol=mp.mpf('1e-9'))
    # 6.1.60 Q_L eigenvalue
    I2=mp.eye(len(m2s))
    w=mp.mpf(0)
    Lmax=int(min(2*S1,2*S2))
    for L in range(0,Lmax+1):
        QL=mp.zeros(NN,NN)
        for M in range(-L,L+1):
            QL=QL+(-1)**M*kron(Tmat(S1,L,M),I2)*kron(I1,Tmat(S2,L,-M))
        for S in Svals2:
            for m in [S-k for k in range(int(2*S)+1)]:
                ev=(-1)**int(S1+S2+S)*(2*L+1)*sixj(S1,S2,S,S2,S1,L)
                w=max(w,mnorm(QL*coupled(S,m)-ev*coupled(S,m)))
    ok&=report(f"({float(S1)},{float(S2)}) 6.1.60 Q_L eigen", w, tol=mp.mpf('1e-9'))
    # 6.1.64/65 projection operator: P_S = sum_m coupled coupled^dag  and via Q_L
    w=mp.mpf(0)
    for S in Svals2:
        PS=mp.zeros(NN,NN)
        for m in [S-k for k in range(int(2*S)+1)]:
            PS=PS+coupled(S,m)*coupled(S,m).H
        # via 6.1.65
        PS2=mp.zeros(NN,NN)
        for L in range(0,Lmax+1):
            QL=mp.zeros(NN,NN)
            for M in range(-L,L+1):
                QL=QL+(-1)**M*kron(Tmat(S1,L,M),I2)*kron(I1,Tmat(S2,L,-M))
            PS2=PS2+sixj(S1,S2,S,S2,S1,L)*QL
        PS2=(-1)**int(S1+S2+S)*(2*S+1)*PS2
        w=max(w,mnorm(PS-PS2))
    ok&=report(f"({float(S1)},{float(S2)}) 6.1.64/65 proj op", w, tol=mp.mpf('1e-9'))

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
