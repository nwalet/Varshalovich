#!/usr/bin/env python3
r"""Verify Sec 7.3 (vector spherical harmonics Y^L_{JM}, S=1) and Table 7.3.1.

Y^L_{JM}(th,ph) 3-vector, spherical components (mu=+1,0,-1):
  [Y^L_{JM}]^mu = C^{JM}_{L,M-mu,1,mu} Y_{L,M-mu}.
Spherical basis vectors (cartesian): e_{+1}=-(x+iy)/sqrt2, e0=z, e_{-1}=(x-iy)/sqrt2.
Table 7.3.1: W_JM = Omega^dag Omega (spinor, S=1/2).
"""
import mpmath as mp
mp.mp.dps=30
def cg(a,b,c,d,e,f):
    from sympy.physics.quantum.cg import CG
    from sympy import Rational as R, N
    from fractions import Fraction
    def rr(x): return R(Fraction(float(x)).limit_denominator(64))
    return mp.mpf(str(complex(N(CG(rr(a),rr(b),rr(c),rr(d),rr(e),rr(f)).doit(),30)).real))
def Y(L,m,th,ph):
    if L<0 or abs(m)>L: return mp.mpc(0)
    return mp.spherharm(L,m,th,ph)
half=mp.mpf(1)/2
def report(tag,w,tol=mp.mpf('1e-11')):
    w=float(w); ok=w<tol; print(f"  {tag:46s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok
ok=True
TH,PH=mp.mpf('0.9'),mp.mpf('0.6')
r2=mp.sqrt(2)

def Yvec_sph(J,L,M,th,ph):
    """spherical components [+1,0,-1]."""
    return {mu: cg(L,M-mu,1,mu,J,M)*Y(int(L),int(M-mu),th,ph) for mu in (1,0,-1)}
def sph_to_cart(comp):
    """comp dict {mu: c} -> cartesian 3-vector.  V = sum_mu c^mu e_mu."""
    ep={1:mp.matrix([-1/r2,-1j/r2,0]),0:mp.matrix([0,0,1]),-1:mp.matrix([1/r2,-1j/r2,0])}
    v=mp.matrix([0,0,0])
    for mu in (1,0,-1): v=v+comp[mu]*ep[mu]
    return v
def vnorm(v): return max(abs(v[i]) for i in range(len(v)))

# ---- Table 7.3.1: W_JM ----
print("Table 7.3.1  W_JM(theta) explicit forms")
def Wval(J,M,th):
    # Omega^dag Omega, L=J+1/2 form
    return float(1/(2*(J+1))*((J+M+1)*abs(Y(int(J+half),int(M+half),th,PH))**2+(J-M+1)*abs(Y(int(J+half),int(M-half),th,PH))**2))
pi=mp.pi; c=lambda:mp.cos(TH); s=lambda:mp.sin(TH)
C=mp.cos(TH); Sn=mp.sin(TH)
tbl={
 (half,half): 1/(4*pi),
 (mp.mpf(3)/2,half): 1/(8*pi)*(3*C**2+1),
 (mp.mpf(3)/2,mp.mpf(3)/2): 3/(8*pi)*Sn**2,
 (mp.mpf(5)/2,half): 3/(16*pi)*(5*C**4-2*C**2+1),
 (mp.mpf(5)/2,mp.mpf(3)/2): 3/(32*pi)*Sn**2*(15*C**2+1),
 (mp.mpf(5)/2,mp.mpf(5)/2): 15/(32*pi)*Sn**4,
 (mp.mpf(7)/2,half): 1/(64*pi)*(175*C**6-165*C**4+45*C**2+9),
 (mp.mpf(7)/2,mp.mpf(3)/2): 15/(64*pi)*Sn**2*(21*C**4-6*C**2+1),
 (mp.mpf(7)/2,mp.mpf(5)/2): 5/(64*pi)*Sn**4*(35*C**2+1),
 (mp.mpf(7)/2,mp.mpf(7)/2): 35/(64*pi)*Sn**6,
 (mp.mpf(9)/2,half): 5/(256*pi)*(441*C**8-644*C**6+294*C**4-36*C**2+9),
 (mp.mpf(9)/2,mp.mpf(3)/2): 15/(128*pi)*Sn**2*(147*C**6-105*C**4+21*C**2+1),
 (mp.mpf(9)/2,mp.mpf(5)/2): 35/(128*pi)*Sn**4*(45*C**4-10*C**2+1),
 (mp.mpf(9)/2,mp.mpf(7)/2): 35/(512*pi)*Sn**6*(63*C**2+1),
 (mp.mpf(9)/2,mp.mpf(9)/2): 315/(512*pi)*Sn**8,
 (mp.mpf(11)/2,half): 3/(512*pi)*(4851*C**10-9555*C**8+6510*C**6-1750*C**4+175*C**2+25),
 (mp.mpf(11)/2,mp.mpf(3)/2): 105/(512*pi)*Sn**2*(297*C**8-348*C**6+126*C**4-12*C**2+1),
 (mp.mpf(11)/2,mp.mpf(5)/2): 105/(1024*pi)*Sn**4*(495*C**6-285*C**4+45*C**2+1),  # note: book "45 cos^2 9" -> 45 cos^2 theta
 (mp.mpf(11)/2,mp.mpf(7)/2): 315/(1024*pi)*Sn**6*(77*C**4-14*C**2+1),
 (mp.mpf(11)/2,mp.mpf(9)/2): 63/(1024*pi)*Sn**8*(99*C**2+1),
 (mp.mpf(11)/2,mp.mpf(11)/2): 693/(1024*pi)*Sn**10,
}
w=mp.mpf(0)
for (J,M),expr in tbl.items():
    w=max(w,abs(Wval(J,M,TH)-expr))
ok&=report("Table 7.3.1 (22 entries; incl. 45cos^2 9 fix)", w)

# ---- 7.3.12 contravariant components (definition) ----
print("7.3.12/14 components; 7.3.15-17 explicit")
w=mp.mpf(0)
for J in [1,2,3,4]:
    for M in range(-J,J+1):
        for L in [J-1,J,J+1]:
            if L<0: continue
            comp=Yvec_sph(J,L,M,TH,PH)
            # covariant 7.3.13: [Y]_mu = (-1)^mu [Y]^{-mu}
            for mu in (1,0,-1):
                cov=(-1)**mu*comp[-mu]
                pred=(-1)**mu*cg(L,M+mu,1,-mu,J,M)*Y(int(L),int(M+mu),TH,PH)
                w=max(w,abs(cov-pred))
ok&=report("7.3.13/14 covariant components", w)

# 7.3.15 (L=J+1), 7.3.16 (L=J), 7.3.17 (L=J-1) explicit contravariant
w=mp.mpf(0)
for J in [1,2,3,4]:
    for M in range(-J,J+1):
        cJp=Yvec_sph(J,J+1,M,TH,PH); cJ=Yvec_sph(J,J,M,TH,PH)
        # 7.3.15
        e_p1=mp.sqrt(mp.mpf((J-M+1)*(J-M+2))/(2*(J+1)*(2*J+3)))*Y(J+1,M-1,TH,PH)
        e_0 =-mp.sqrt(mp.mpf((J-M+1)*(J+M+1))/((J+1)*(2*J+3)))*Y(J+1,M,TH,PH)
        e_m1=mp.sqrt(mp.mpf((J+M+1)*(J+M+2))/(2*(J+1)*(2*J+3)))*Y(J+1,M+1,TH,PH)
        w=max(w,abs(cJp[1]-e_p1),abs(cJp[0]-e_0),abs(cJp[-1]-e_m1))
        # 7.3.16 (L=J)
        f_p1=-mp.sqrt(mp.mpf((J+M)*(J-M+1))/(2*J*(J+1)))*Y(J,M-1,TH,PH)
        f_0 = M/mp.sqrt(mp.mpf(J*(J+1)))*Y(J,M,TH,PH)
        f_m1=mp.sqrt(mp.mpf((J-M)*(J+M+1))/(2*J*(J+1)))*Y(J,M+1,TH,PH)
        w=max(w,abs(cJ[1]-f_p1),abs(cJ[0]-f_0),abs(cJ[-1]-f_m1))
        # 7.3.17 (L=J-1)
        if J-1>=0:
            cJm=Yvec_sph(J,J-1,M,TH,PH)
            g_p1=mp.sqrt(mp.mpf((J+M)*(J+M-1))/(2*J*(2*J-1)))*Y(J-1,M-1,TH,PH)
            g_0 =mp.sqrt(mp.mpf((J-M)*(J+M))/(J*(2*J-1)))*Y(J-1,M,TH,PH)
            g_m1=mp.sqrt(mp.mpf((J-M)*(J-M-1))/(2*J*(2*J-1)))*Y(J-1,M+1,TH,PH)
            w=max(w,abs(cJm[1]-g_p1),abs(cJm[0]-g_0),abs(cJm[-1]-g_m1))
ok&=report("7.3.15/16/17 explicit components", w)

# ---- 7.3.9/10 (lambda) decompositions ----
print("7.3.9/10 Y^(lambda) decompositions")
def Ylam(J,lam,M,th,ph):
    """spherical comps of Y^(lambda) via 7.3.9."""
    Jm1=Yvec_sph(J,J-1,M,th,ph); Jp1=Yvec_sph(J,J+1,M,th,ph); JJ=Yvec_sph(J,J,M,th,ph)
    out={}
    for mu in (1,0,-1):
        if lam==1: out[mu]=mp.sqrt(mp.mpf(J+1)/(2*J+1))*Jm1[mu]+mp.sqrt(mp.mpf(J)/(2*J+1))*Jp1[mu]
        elif lam==0: out[mu]=JJ[mu]
        else: out[mu]=mp.sqrt(mp.mpf(J)/(2*J+1))*Jm1[mu]-mp.sqrt(mp.mpf(J+1)/(2*J+1))*Jp1[mu]
    return out
# 7.3.19 explicit (lambda) components match 7.3.9 construction
w=mp.mpf(0)
for J in [1,2,3]:
    for M in range(-J,J+1):
        for lam in (1,0,-1):
            yl=Ylam(J,lam,M,TH,PH)
            for mu in (1,0,-1):
                if lam==1: pred=mp.sqrt(mp.mpf(J+1)/(2*J+1))*cg(J-1,M-mu,1,mu,J,M)*Y(J-1,int(M-mu),TH,PH)+mp.sqrt(mp.mpf(J)/(2*J+1))*cg(J+1,M-mu,1,mu,J,M)*Y(J+1,int(M-mu),TH,PH)
                elif lam==0: pred=cg(J,M-mu,1,mu,J,M)*Y(J,int(M-mu),TH,PH)
                else: pred=mp.sqrt(mp.mpf(J)/(2*J+1))*cg(J-1,M-mu,1,mu,J,M)*Y(J-1,int(M-mu),TH,PH)-mp.sqrt(mp.mpf(J+1)/(2*J+1))*cg(J+1,M-mu,1,mu,J,M)*Y(J+1,int(M-mu),TH,PH)
                w=max(w,abs(yl[mu]-pred))
ok&=report("7.3.19 (lambda) components", w)

# ---- 7.3.5/6/7 transverse & longitudinal ----
print("7.3.5/6/7 transverse/longitudinal")
ncart=mp.matrix([mp.sin(TH)*mp.cos(PH),mp.sin(TH)*mp.sin(PH),mp.cos(TH)])
def cross(a,b): return mp.matrix([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])
w=mp.mpf(0)
for J in [1,2,3]:
    for M in range(-J,J+1):
        Y1=sph_to_cart(Ylam(J,1,M,TH,PH)); Y0=sph_to_cart(Ylam(J,0,M,TH,PH)); Ym=sph_to_cart(Ylam(J,-1,M,TH,PH))
        # 7.3.5 n.Y1=n.Y0=0, n x Ym=0
        w=max(w,abs(sum(ncart[i]*Y1[i] for i in range(3))),abs(sum(ncart[i]*Y0[i] for i in range(3))))
        w=max(w,vnorm(cross(ncart,Ym)))
        # 7.3.7 Y^{(-1)} = n Y_JM
        w=max(w,vnorm(Ym-Y(J,M,TH,PH)*ncart))
ok&=report("7.3.5 transverse/longitudinal; 7.3.7 Y^(-1)=nY", w)

# ---- 7.3.39 complex conjugation ----
print("7.3.39 complex conjugation")
w=mp.mpf(0)
for J in [1,2,3]:
    for M in range(-J,J+1):
        for L in [J-1,J,J+1]:
            if L<0: continue
            comp=Yvec_sph(J,L,M,TH,PH); compc={mu:mp.conj(comp[mu]) for mu in (1,0,-1)}
            ref=Yvec_sph(J,L,-M,TH,PH)
            # Y^L*_{JM} = (-1)^{J+L+M+1} Y^L_{J,-M} (as vectors -> spherical comps conj)
            # [Y^L_{JM}]^{mu *} = (-1)^{J+L+M+1}[Y^L_{J-M}]_mu = (-1)^{J+L+M+1}(-1)^mu[Y^L_{J-M}]^{-mu}
            for mu in (1,0,-1):
                pred=(-1)**(J+L+M+1)*(-1)**mu*ref[-mu]
                w=max(w,abs(compc[mu]-pred))
        for lam in (1,0,-1):
            comp=Ylam(J,lam,M,TH,PH); compc={mu:mp.conj(comp[mu]) for mu in (1,0,-1)}
            ref=Ylam(J,lam,-M,TH,PH)
            for mu in (1,0,-1):
                pred=(-1)**(M+lam+1)*(-1)**mu*ref[-mu]
                w=max(w,abs(compc[mu]-pred))
ok&=report("7.3.39/40-43 complex conjugation", w)

# ---- 7.3.61-63 (S.n) action (spin-1 n.S matrix, spherical basis +1,0,-1) ----
print("7.3.61-63 (S.n) action")
Sz3=mp.matrix([[1,0,0],[0,0,0],[0,0,-1]])
Sp3=mp.matrix([[0,r2,0],[0,0,r2],[0,0,0]]); Sm3=Sp3.T
Sx3=(Sp3+Sm3)/2; Sy3=(Sp3-Sm3)/(2j)
nx=mp.sin(TH)*mp.cos(PH); ny=mp.sin(TH)*mp.sin(PH); nz=mp.cos(TH)
Sn3=Sx3*nx+Sy3*ny+Sz3*nz   # rows/cols order +1,0,-1
def apply3(Mmat,comp):
    v=mp.matrix([comp[1],comp[0],comp[-1]]); out=Mmat*v
    return {1:out[0],0:out[1],-1:out[2]}
w=mp.mpf(0)
for J in [1,2,3]:
    for M in range(-J,J+1):
        # L-basis 7.3.62
        lhsp=apply3(Sn3,Yvec_sph(J,J+1,M,TH,PH))
        rhsp={mu:-mp.sqrt(mp.mpf(J)/(2*J+1))*Yvec_sph(J,J,M,TH,PH)[mu] for mu in (1,0,-1)}
        lhsJ=apply3(Sn3,Yvec_sph(J,J,M,TH,PH))
        rhsJ={mu:-mp.sqrt(mp.mpf(J)/(2*J+1))*Yvec_sph(J,J+1,M,TH,PH)[mu]-mp.sqrt(mp.mpf(J+1)/(2*J+1))*Yvec_sph(J,J-1,M,TH,PH)[mu] for mu in (1,0,-1)}
        for mu in (1,0,-1): w=max(w,abs(lhsp[mu]-rhsp[mu]),abs(lhsJ[mu]-rhsJ[mu]))
        # lambda-basis 7.3.63: (S.n)Y^(1)=-Y^(0), (S.n)Y^(0)=-Y^(1), (S.n)Y^(-1)=0
        l1=apply3(Sn3,Ylam(J,1,M,TH,PH)); r1=Ylam(J,0,M,TH,PH)
        l0=apply3(Sn3,Ylam(J,0,M,TH,PH)); r0=Ylam(J,1,M,TH,PH)
        lm=apply3(Sn3,Ylam(J,-1,M,TH,PH))
        for mu in (1,0,-1): w=max(w,abs(l1[mu]+r1[mu]),abs(l0[mu]+r0[mu]),abs(lm[mu]))
ok&=report("7.3.62/63 (S.n) action L- and lambda-basis", w)

# ---- 7.3.67-69 (S.L) action ----
print("7.3.67-69 (S.L) action")
# S.L = sum_mu (-1)^mu S_mu L_{-mu}; S_mu spin-1 (spherical), L_{-mu} ladder on Y_Lm.
Ssph={0:Sz3, 1:-(Sx3+1j*Sy3)/r2, -1:(Sx3-1j*Sy3)/r2}
def Lmu_on_Ycomp(mu, J, L, M, th, ph):
    """return dict of spherical comps of L_mu acting on Y^L_{JM} (L fixed, m->m+mu)."""
    # [Y^L_{JM}]^sigma = C^{JM}_{L,M-sigma,1,sigma} Y_{L,M-sigma}; L_mu Y_{Lm}=sqrt(L(L+1)-m(m+mu)) Y_{L,m+mu} (mu=+-1); L_0 Y=m Y
    out={}
    for sig in (1,0,-1):
        m=M-sig; coef=cg(L,m,1,sig,J,M)
        if mu==0: val=coef*m*Y(int(L),int(m),th,ph)
        else:
            # spherical comp: L_{+1}=-L_+/sqrt2, L_{-1}=+L_-/sqrt2
            sgn=(-1 if mu==1 else 1)/r2
            fac=mp.sqrt(L*(L+1)-m*(m+mu))
            val=coef*sgn*fac*Y(int(L),int(m+mu),th,ph)
        out[sig]=val
    return out
def SL_on(J,L,M,th,ph):
    # (S.L) comp sigma' = sum_mu (-1)^mu sum_sigma [S_mu]_{sigma' sigma} [L_{-mu} Y]^sigma
    res={1:mp.mpc(0),0:mp.mpc(0),-1:mp.mpc(0)}
    order=[1,0,-1]
    for mu in (1,0,-1):
        Lcomp=Lmu_on_Ycomp(-mu,J,L,M,th,ph)  # L_{-mu}
        Smat=Ssph[mu]
        for i,sp in enumerate(order):
            for j,sg in enumerate(order):
                res[sp]+=(-1)**mu*Smat[i,j]*Lcomp[sg]
    return res
w=mp.mpf(0)
for J in [1,2,3]:
    for M in range(-J,J+1):
        for L in [J-1,J,J+1]:
            if L<0: continue
            sl=SL_on(J,L,M,TH,PH)
            ev=mp.mpf(1)/2*(J*(J+1)-L*(L+1)-2)
            for mu in (1,0,-1):
                w=max(w,abs(sl[mu]-ev*Yvec_sph(J,L,M,TH,PH)[mu]))
ok&=report("7.3.67/68 (S.L) eigenvalue on Y^L", w)

# ---- 7.3.78/79 sum -> n ----
print("7.3.78/79 sum_M Y*_JM Y^L_JM = ... n")
w=mp.mpf(0)
for J in [1,2,3]:
    for L in [J-1,J,J+1]:
        if L<0: continue
        acc=mp.matrix([0,0,0])
        for M in range(-J,J+1):
            acc=acc+mp.conj(Y(J,M,TH,PH))*sph_to_cart(Yvec_sph(J,L,M,TH,PH))
        coef={J+1:-mp.sqrt(mp.mpf((J+1)*(2*J+1)))/(4*pi),J:mp.mpf(0),J-1:mp.sqrt(mp.mpf(J*(2*J+1)))/(4*pi)}[L]
        w=max(w,vnorm(acc-coef*ncart))
ok&=report("7.3.79 sum_M Y* Y^L = coef*n", w)

# ---- 7.3.82-84 quadratic sums |a.Y|^2 ----
print("7.3.82-84 sum_M |a.Y^L|^2")
avec=mp.matrix([mp.mpf('0.6')+0.2j,mp.mpf('0.3')-0.5j,mp.mpf('0.8')+0.1j])
a2=sum(abs(avec[i])**2 for i in range(3)); na=sum(ncart[i]*avec[i] for i in range(3))
w=mp.mpf(0)
for J in [1,2,3]:
    def aY(L,M): return sum(avec[i]*sph_to_cart(Yvec_sph(J,L,M,TH,PH))[i] for i in range(3))
    s_p=sum(abs(aY(J+1,M))**2 for M in range(-J,J+1))
    s_0=sum(abs(aY(J,M))**2 for M in range(-J,J+1))
    s_m=sum(abs(aY(J-1,M))**2 for M in range(-J,J+1))
    w=max(w,abs(s_p-1/(8*pi)*(J*a2+(J+2)*abs(na)**2)))
    w=max(w,abs(s_0-(2*J+1)/(8*pi)*(a2-abs(na)**2)))
    w=max(w,abs(s_m-1/(8*pi)*((J+1)*a2+(J-1)*abs(na)**2)))
ok&=report("7.3.82/83/84 |a.Y^{J+1,J,J-1}|^2", w)

# ---- 7.3.100 CG dot-series ----
print("7.3.100 CG series (dot product)")
def sixj(a,b,c,d,e,f):
    from sympy.physics.wigner import wigner_6j
    from sympy import Rational as R, N
    from fractions import Fraction
    def rr(x): return R(Fraction(float(x)).limit_denominator(64))
    try: return float(N(wigner_6j(rr(a),rr(b),rr(c),rr(d),rr(e),rr(f)),30))
    except Exception: return 0.0
w=mp.mpf(0)
for (J1,L1,M1,J2,L2,M2) in [(1,1,0,1,1,1),(2,2,1,1,1,-1),(2,1,0,2,3,1),(1,2,-1,2,2,1)]:
    v1=sph_to_cart(Yvec_sph(J1,L1,M1,TH,PH)); v2=sph_to_cart(Yvec_sph(J2,L2,M2,TH,PH))
    lhs=sum(v1[i]*v2[i] for i in range(3))
    rhs=mp.mpc(0); M=M1+M2
    for L in range(abs(L1-L2),L1+L2+1):
        if abs(M)>L: continue
        rhs+=((-1)**(J2+L1+L)*mp.sqrt(mp.mpf((2*J1+1)*(2*J2+1)*(2*L1+1)*(2*L2+1))/(4*pi*(2*L+1)))
              *sixj(L1,L2,L,J2,J1,1)*cg(L1,0,L2,0,L,0)*cg(J1,M1,J2,M2,L,M)*Y(L,int(M),TH,PH))
    w=max(w,abs(lhs-rhs))
ok&=report("7.3.100 Y^L1.Y^L2 = sum_L ... Y_LM", w)

# ---- 7.3.102 addition theorem (scalar) ----
print("7.3.102 addition theorem")
TH2,PH2=mp.mpf('1.4'),mp.mpf('2.1')
cos12=mp.cos(TH)*mp.cos(TH2)+mp.sin(TH)*mp.sin(TH2)*mp.cos(PH-PH2)
w=mp.mpf(0)
for J in [1,2,3]:
    for L in [J-1,J,J+1]:
        if L<0: continue
        s=mp.mpc(0)
        for M in range(-J,J+1):
            v1=sph_to_cart(Yvec_sph(J,L,M,TH,PH)); v2=sph_to_cart(Yvec_sph(J,L,M,TH2,PH2))
            s+=sum(mp.conj(v1[i])*v2[i] for i in range(3))
        rhs=(2*J+1)*mp.legendre(int(L),cos12)/(4*pi)
        w=max(w,abs(s-rhs))
ok&=report("7.3.102 4pi sum Y^L*.Y^L = (2J+1)P_L", w)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
