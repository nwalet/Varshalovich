#!/usr/bin/env python3
r"""Verify Sec 7.2 (spinor spherical harmonics Omega^L_{JM}, S=1/2).

Omega^L_{JM}(th,ph) is a 2-spinor, contravariant components
  [Omega]^mu = C^{J M}_{L, M-mu, 1/2, mu} Y_{L,M-mu}(th,ph),  mu=+1/2,-1/2.
Basis order (row 0,1) = (mu=+1/2, mu=-1/2).  Y = mpmath.spherharm (VMK conv).
clebsch{a}{b}{c}{d}{e}{f}=C^{e f}_{a b,c d}.
"""
import mpmath as mp
mp.mp.dps = 30
def cg(a,b,c,d,e,f):
    from sympy.physics.quantum.cg import CG
    from sympy import Rational as R, N
    from fractions import Fraction
    def rr(x): return R(Fraction(float(x)).limit_denominator(64))
    return mp.mpf(str(complex(N(CG(rr(a),rr(b),rr(c),rr(d),rr(e),rr(f)).doit(),30)).real))
def Y(L,m,th,ph):
    if abs(m)>L: return mp.mpc(0)
    return mp.spherharm(L,m,th,ph)
half=mp.mpf(1)/2
def Omega(J,L,M,th,ph):
    """2-vector [Om^{+1/2}, Om^{-1/2}]."""
    out=[]
    for mu in (half,-half):
        out.append(cg(L,M-mu,half,mu,J,M)*Y(L,int(M-mu),th,ph))
    return mp.matrix(out)
# spin-1/2
sx=mp.matrix([[0,1],[1,0]])/2; sy=mp.matrix([[0,-1j],[1j,0]])/2; sz=mp.matrix([[1,0],[0,-1]])/2
def report(tag,w,tol=mp.mpf('1e-12')):
    w=float(w); ok=w<tol; print(f"  {tag:48s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok
def vnorm(v): return max(abs(v[i]) for i in range(len(v)))
ok=True
TH,PH=mp.mpf('1.0'),mp.mpf('0.7')
# (J,L) pairs with L=J+-1/2
JLs=[(half,1),(half,0),(mp.mpf(3)/2,2),(mp.mpf(3)/2,1),(mp.mpf(5)/2,3),(mp.mpf(5)/2,2),(mp.mpf(7)/2,4),(mp.mpf(7)/2,3)]

# ---- 7.2.6/7 explicit contravariant forms ----
print("7.2.6/7 explicit L=J+-1/2 contravariant components")
w=mp.mpf(0)
for J in [half,mp.mpf(3)/2,mp.mpf(5)/2,mp.mpf(7)/2]:
    for M in [J-k for k in range(int(2*J)+1)]:
        # L=J+1/2
        Lp=J+half; om=Omega(J,Lp,M,TH,PH)
        e1=-mp.sqrt((J-M+1)/(2*(J+1)))*Y(int(Lp),int(M-half),TH,PH)
        e2= mp.sqrt((J+M+1)/(2*(J+1)))*Y(int(Lp),int(M+half),TH,PH)
        w=max(w,abs(om[0]-e1),abs(om[1]-e2))
        # L=J-1/2
        Lm=J-half
        if Lm>=0:
            om=Omega(J,Lm,M,TH,PH)
            e1=mp.sqrt((J+M)/(2*J))*Y(int(Lm),int(M-half),TH,PH)
            e2=mp.sqrt((J-M)/(2*J))*Y(int(Lm),int(M+half),TH,PH)
            w=max(w,abs(om[0]-e1),abs(om[1]-e2))
ok&=report("7.2.6/7 explicit forms", w)

# ---- 7.2.9/10/11 covariant components [Om]_mu = (-1)^{1/2-mu} [Om]^{-mu} ----
print("7.2.9-11 covariant components")
w=mp.mpf(0)
for J in [half,mp.mpf(3)/2,mp.mpf(5)/2]:
    for M in [J-k for k in range(int(2*J)+1)]:
        for L in ([J+half]+([J-half] if J-half>=0 else [])):
            om=Omega(J,L,M,TH,PH)
            cov={}
            for i,mu in enumerate((half,-half)):
                cov[mu]=(-1)**int(half-mu)*om[1 if mu==half else 0]  # [Om]^{-mu}
            # 7.2.9 formula
            for mu in (half,-half):
                pred=(-1)**int(half-mu)*cg(L,M+mu,half,-mu,J,M)*Y(int(L),int(M+mu),TH,PH)
                w=max(w,abs(cov[mu]-pred))
ok&=report("7.2.9 covariant = (-1)^{1/2-mu} C Y", w)

# ---- 7.2.14 helicity expansion ----
print("7.2.14 helicity expansion vs direct")
# need helicity spin functions chi_{1/2 lam}(th,ph) and D^{1/2}? No: 7.2.14 uses D^J.
def d1half(b):
    c=mp.cos(b/2); s=mp.sin(b/2); return mp.matrix([[c,-s],[s,c]])
def helchi(lam,th,ph):  # chi_{1/2 lam}(th,ph) column, from 6.2.22
    if lam==half: return mp.matrix([mp.cos(th/2)*mp.e**(-1j*ph/2), mp.sin(th/2)*mp.e**(1j*ph/2)])
    else: return mp.matrix([-mp.sin(th/2)*mp.e**(-1j*ph/2), mp.cos(th/2)*mp.e**(1j*ph/2)])
def DJ(J,Mrow,Mcol,a,b,g):
    # D^J_{Mrow,Mcol}(a,b,g)=e^{-i Mrow a} d^J e^{-i Mcol g}; build d^J via wigner_d
    import sys,os; sys.path.insert(0,os.path.dirname(__file__))
    from wigner_d import wigner_d, beta as B
    d=complex(wigner_d(float(J),float(Mrow),float(Mcol)).subs(B,mp.mpf(b)).evalf(25))
    return mp.e**(-1j*Mrow*a)*mp.mpf(d.real)*mp.e**(-1j*Mcol*g)
w=mp.mpf(0)
for J in [half,mp.mpf(3)/2,mp.mpf(5)/2]:
    for M in [J-k for k in range(int(2*J)+1)]:
        for L in ([J+half]+([J-half] if J-half>=0 else [])):
            direct=Omega(J,L,M,TH,PH)
            hel=mp.matrix([0,0])
            for lam in (half,-half):
                hel=hel+mp.sqrt((2*L+1)/(4*mp.pi))*cg(L,0,half,lam,J,lam)*DJ(J,-lam,-M,0,TH,PH)*helchi(lam,TH,PH)
            w=max(w,vnorm(direct-hel))
ok&=report("7.2.14 helicity expansion", w)

# ---- 7.2.16 complex conjugation: Om* = (-1)^{J+L-M} i sigma_y Om_{J,-M} ----
print("7.2.16 complex conjugation")
isy=mp.matrix([[0,1],[-1,0]])
w=mp.mpf(0)
for J in [half,mp.mpf(3)/2,mp.mpf(5)/2]:
    for M in [J-k for k in range(int(2*J)+1)]:
        for L in ([J+half]+([J-half] if J-half>=0 else [])):
            om=Omega(J,L,M,TH,PH); omc=mp.matrix([mp.conj(om[i]) for i in range(2)])
            rhs=(-1)**int(J+L-M)*(isy*Omega(J,L,-M,TH,PH))
            w=max(w,vnorm(omc-rhs))
ok&=report("7.2.16 Om* = (-1)^{J+L-M} i sy Om_{-M}", w)

# ---- 7.2.20 time reversal: sy Om* = (-1)^{J+L-M+1/2} Om_{-M} ----
print("7.2.20 time reversal")
sy_m=mp.matrix([[0,-1j],[1j,0]])
w=mp.mpf(0)
for J in [half,mp.mpf(3)/2,mp.mpf(5)/2]:
    for M in [J-k for k in range(int(2*J)+1)]:
        for L in ([J+half]+([J-half] if J-half>=0 else [])):
            om=Omega(J,L,M,TH,PH); omc=mp.matrix([mp.conj(om[i]) for i in range(2)])
            lhs=sy_m*omc
            rhs=mp.e**(1j*mp.pi*(J+L-M+half))*Omega(J,L,-M,TH,PH)
            w=max(w,vnorm(lhs-rhs))
ok&=report("7.2.20 sy Om* = (-1)^{J+L-M+1/2} Om_{-M}", w)

# ---- 7.2.23/25 (S.n) Om^L = -1/2 Om^{L'}, L'=2J-L ----
print("7.2.23/25 (S.n) action")
nx=mp.sin(TH)*mp.cos(PH); ny=mp.sin(TH)*mp.sin(PH); nz=mp.cos(TH)
Sn=sx*nx+sy*ny+sz*nz
w=mp.mpf(0)
for J in [half,mp.mpf(3)/2,mp.mpf(5)/2]:
    for M in [J-k for k in range(int(2*J)+1)]:
        for L in ([J+half]+([J-half] if J-half>=0 else [])):
            Lp=int(2*J-L)
            lhs=Sn*Omega(J,L,M,TH,PH)
            rhs=-half*Omega(J,Lp,M,TH,PH)
            w=max(w,vnorm(lhs-rhs))
ok&=report("7.2.23/25 (S.n) Om = -1/2 Om^{2J-L}", w)

# ---- 7.2.29/30 (L.S) eigenvalue = 1/2{J(J+1)-L(L+1)-3/4} ----
print("7.2.29/30 (L.S) eigenvalue (via S^2 identity)")
# (L.S) = 1/2(J^2 - L^2 - S^2); J,L,S are quantum numbers -> pure eigenvalue check
w=mp.mpf(0)
for J in [half,mp.mpf(3)/2,mp.mpf(5)/2,mp.mpf(7)/2]:
    for L in ([J+half]+([J-half] if J-half>=0 else [])):
        ev=half*(J*(J+1)-L*(L+1)-mp.mpf(3)/4)
        # detailed 7.2.30: J+1/2 -> -(2J+3)/4 ; J-1/2 -> (2J-1)/4
        if L==J+half: pred=-(2*J+3)/4
        else: pred=(2*J-1)/4
        w=max(w,abs(ev-pred))
ok&=report("7.2.30 L.S eigenvalues", w)

# ---- 7.2.46/47 quadratic form W_JM = Om^dag Om ----
print("7.2.46/47 W_JM = Om^dag Om (L-independent)")
w=mp.mpf(0)
for J in [half,mp.mpf(3)/2,mp.mpf(5)/2]:
    for M in [J-k for k in range(int(2*J)+1)]:
        oms=[]
        for L in ([J+half]+([J-half] if J-half>=0 else [])):
            om=Omega(J,L,M,TH,PH)
            oms.append(sum(abs(om[i])**2 for i in range(2)))
        # W from 7.2.47 (first form)
        Wpred=1/(2*(J+1))*((J+M+1)*abs(Y(int(J+half),int(M+half),TH,PH))**2+(J-M+1)*abs(Y(int(J+half),int(M-half),TH,PH))**2)
        for wv in oms: w=max(w,abs(wv-Wpred))
        # second form (L=J-1/2)
        if J-half>=0:
            W2=1/(2*J)*((J+M)*abs(Y(int(J-half),int(M-half),TH,PH))**2+(J-M)*abs(Y(int(J-half),int(M+half),TH,PH))**2)
            w=max(w,abs(W2-Wpred))
ok&=report("7.2.46/47 W_JM both forms, L-independent", w)

# ---- 7.2.50 low Legendre coeffs a0,a1,a2 via integration ----
print("7.2.50 Legendre coefficients a0,a1,a2")
def Wfun(J,M,th):
    return float(1/(2*(J+1))*((J+M+1)*abs(Y(int(J+half),int(M+half),th,mp.mpf('0.3')))**2+(J-M+1)*abs(Y(int(J+half),int(M-half),th,mp.mpf('0.3')))**2))
def an(J,M,n):
    f=lambda th: Wfun(J,M,th)*float(mp.legendre(2*n,mp.cos(th)))*mp.sin(th)
    integ=mp.quad(f,[0,mp.pi])
    return (4*n+1)/2*integ   # a_m = (4m+1)/2 int_0^pi W P_{2m}(cos th) sin th dth
w=mp.mpf(0)
for J in [mp.mpf(3)/2,mp.mpf(5)/2,mp.mpf(7)/2]:
    for M in [J-k for k in range(int(2*J)+1) if (J-k)>0]:
        a0=an(J,M,0); a1=an(J,M,1)
        p0=1/(4*mp.pi); p1=mp.mpf(5)/(16*mp.pi)*(J*(J+1)-3*M*M)/(J*(J+1))
        w=max(w,abs(a0-p0),abs(a1-p1))
ok&=report("7.2.50 a0=1/4pi, a1", w, tol=mp.mpf('1e-9'))

# ---- 7.2.31 cos-theta recursion ----
print("7.2.31 cos(theta) recursion")
def Om(J,L,M):
    if J<half or L<0 or abs(M)>J or L not in (J+half,J-half): return mp.matrix([0,0])
    return Omega(J,L,M,TH,PH)
w=mp.mpf(0)
for J in [mp.mpf(3)/2,mp.mpf(5)/2]:
    for M in [J-k for k in range(int(2*J)+1)]:
        for L in ([J+half]+([J-half] if J-half>=0 else [])):
            Lp=2*J-L
            lhs=mp.cos(TH)*Om(J,L,M)
            rhs=(mp.sqrt((J-M+1)*(J+M+1))/(2*(J+1))*Om(J+1,L+1,M)
                 -M/(2*J*(J+1))*Om(J,Lp,M)
                 +mp.sqrt((J-M)*(J+M))/(2*J)*Om(J-1,L-1,M))
            w=max(w,vnorm(lhs-rhs))
ok&=report("7.2.31 cos-theta recursion", w)

# ---- 7.2.40 CG series: verify phase exponent 1/2 (flagged anomaly, was 1/3) ----
print("7.2.40 CG series (Om^dag Om = sum_L ... Y_LM); phase 1/2 vs 1/3")
def sixj(a,b,c,d,e,f):
    from sympy.physics.wigner import wigner_6j
    from sympy import Rational as R, N
    from fractions import Fraction
    def rr(x): return R(Fraction(float(x)).limit_denominator(64))
    try: return float(N(wigner_6j(rr(a),rr(b),rr(c),rr(d),rr(e),rr(f)),30))
    except Exception: return 0.0
def cgseries(J1,L1,M1,J2,L2,M2,phase_half):
    o1=Omega(J1,L1,M1,TH,PH); o2=Omega(J2,L2,M2,TH,PH)
    lhs=sum(mp.conj(o1[i])*o2[i] for i in range(2))
    ph=half if phase_half else mp.mpf(1)/3
    rhs=mp.mpc(0)
    for L in range(abs(L1-L2),L1+L2+1):
        M=M2-M1
        if abs(M)>L: continue
        rhs+=((-1)**(J1+M1+J2+L+ph if False else 0))  # placeholder; compute phase below
    # compute properly (phase may be half-integer -> use e^{i pi x})
    rhs=mp.mpc(0)
    for L in range(abs(L1-L2),L1+L2+1):
        M=M2-M1
        if abs(M)>L: continue
        phase=mp.e**(1j*mp.pi*(J1+M1+J2+L+ph))
        rhs+=(phase*sixj(L1,L2,L,J2,J1,half)
              *mp.sqrt((2*J1+1)*(2*J2+1)*(2*L1+1)*(2*L2+1)/(4*mp.pi*(2*L+1)))
              *cg(L1,0,L2,0,L,0)*cg(J1,-M1,J2,M2,L,M)*Y(L,int(M),TH,PH))
    return abs(lhs-rhs)
w2=mp.mpf(0); w3=mp.mpf(0)
for (J1,L1,J2,L2,M1,M2) in [(half,1,half,1,half,-half),(mp.mpf(3)/2,2,mp.mpf(3)/2,1,half,mp.mpf(3)/2),
                            (mp.mpf(3)/2,1,half,1,-half,half),(mp.mpf(5)/2,3,mp.mpf(3)/2,2,half,-half)]:
    w2=max(w2,cgseries(J1,L1,M1,J2,L2,M2,True))
    w3=max(w3,cgseries(J1,L1,M1,J2,L2,M2,False))
print(f"    phase 1/2: worst={float(w2):.2e}   phase 1/3: worst={float(w3):.2e}")
ok&=report("7.2.40 CG series (phase 1/2 correct)", w2)

# ---- 7.2.42 addition theorem ----
print("7.2.42 addition theorem")
TH2,PH2=mp.mpf('1.4'),mp.mpf('2.1')
cos12=mp.cos(TH)*mp.cos(TH2)+mp.sin(TH)*mp.sin(TH2)*mp.cos(PH-PH2)
w=mp.mpf(0)
for J in [half,mp.mpf(3)/2,mp.mpf(5)/2]:
    for L in ([J+half]+([J-half] if J-half>=0 else [])):
        s=mp.mpc(0)
        for M in [J-k for k in range(int(2*J)+1)]:
            o1=Omega(J,L,M,TH,PH); o2=Omega(J,L,M,TH2,PH2)
            s+=sum(mp.conj(o1[i])*o2[i] for i in range(2))
        rhs=(2*J+1)*mp.legendre(int(L),cos12)/(4*mp.pi)
        w=max(w,abs(s-rhs))
ok&=report("7.2.42 4pi sum Om^dag Om = (2J+1)P_L(cos w12)", w)

# ---- 7.2.54/56 special-M W forms ----
print("7.2.54/56 special-M quadratic forms")
def Wval(J,M):
    return float(1/(2*(J+1))*((J+M+1)*abs(Y(int(J+half),int(M+half),TH,PH))**2+(J-M+1)*abs(Y(int(J+half),int(M-half),TH,PH))**2))
w=mp.mpf(0)
for J in [mp.mpf(3)/2,mp.mpf(5)/2,mp.mpf(7)/2]:
    # 7.2.54 M=+-1/2
    Wm=Wval(J,half)
    Lh=int(J-half)
    Pd=lambda x: float(mp.diff(lambda t: mp.legendre(Lh,t), x))
    pred=1/(2*mp.pi*(2*J+1))*(mp.sin(TH)**2*Pd(mp.cos(TH))**2+(J+half)**2*mp.legendre(Lh,mp.cos(TH))**2)
    w=max(w,abs(Wm-pred))
    # 7.2.56 M=+-J
    WJ=Wval(J,J)
    predJ=mp.factorial(2*J)/(mp.pi*2**(2*J+1)*mp.factorial(J-half)**2)*mp.sin(TH)**(2*J-1)
    w=max(w,abs(WJ-predJ))
ok&=report("7.2.54/56 W_{J,1/2} and W_{J,J}", w, tol=mp.mpf('1e-9'))

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
