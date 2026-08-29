#!/usr/bin/env python3
r"""Verify Sec 7.3.15-17: J=0,1 forms, quadratic forms W^perp/W^par, Table 7.3.2.

W^perp_JM = |Y^(1)_JM|^2 = |Y^(0)_JM|^2 ; W^par_JM = |Y^(-1)_JM|^2 = |Y_JM|^2.
Table entries encoded with OCR corrections; checker flags any wrong reading.
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
r2=mp.sqrt(2); pi=mp.pi
def report(tag,w,tol=mp.mpf('1e-11')):
    w=float(w); ok=w<tol; print(f"  {tag:40s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok
ok=True
def Yvec_sph(J,L,M,th,ph): return {mu:cg(L,M-mu,1,mu,J,M)*Y(int(L),int(M-mu),th,ph) for mu in (1,0,-1)}
def Ylam(J,lam,M,th,ph):
    Jm1=Yvec_sph(J,J-1,M,th,ph); Jp1=Yvec_sph(J,J+1,M,th,ph); JJ=Yvec_sph(J,J,M,th,ph)
    out={}
    for mu in (1,0,-1):
        if lam==1: out[mu]=mp.sqrt(mp.mpf(J+1)/(2*J+1))*Jm1[mu]+mp.sqrt(mp.mpf(J)/(2*J+1))*Jp1[mu]
        elif lam==0: out[mu]=JJ[mu]
        else: out[mu]=mp.sqrt(mp.mpf(J)/(2*J+1))*Jm1[mu]-mp.sqrt(mp.mpf(J+1)/(2*J+1))*Jp1[mu]
    return out
def Wperp(J,M,th): return float(sum(abs(Ylam(J,1,M,th,mp.mpf('0.4'))[mu])**2 for mu in (1,0,-1)))
def Wpar(J,M,th): return float(abs(Y(J,M,th,mp.mpf('0.4')))**2)

TH=mp.mpf('0.8'); C=mp.cos(TH); Sn=mp.sin(TH)

# ---- 7.3.155 definition ----
print("7.3.155 W^perp definition")
w=mp.mpf(0)
for J in [1,2,3,4]:
    for M in range(-J,J+1):
        pred=1/(2*J*(J+1))*((J-M)*(J+M+1)*abs(Y(J,M+1,TH,mp.mpf('0.4')))**2+2*M*M*abs(Y(J,M,TH,mp.mpf('0.4')))**2+(J+M)*(J-M+1)*abs(Y(J,M-1,TH,mp.mpf('0.4')))**2)
        w=max(w,abs(Wperp(J,M,TH)-pred))
        # W^(0) = W^(1)
        w0=sum(abs(Ylam(J,0,M,TH,mp.mpf('0.4'))[mu])**2 for mu in (1,0,-1))
        w=max(w,abs(float(w0)-Wperp(J,M,TH)))
ok&=report("7.3.155 W^perp; W^(0)=W^(1)", w)

# ---- 7.3.159-163 Legendre coeffs a_n,b_n (via integration) ----
print("7.3.159-163 a0,a1,b1 coefficients")
def coeff(Wf,J,M,n):
    f=lambda th: Wf(J,M,th)*float(mp.legendre(2*n,mp.cos(th)))*mp.sin(th)
    return (4*n+1)/2*mp.quad(f,[0,pi])
w=mp.mpf(0)
for J in [2,3,4]:
    for M in range(0,J+1):
        a0=coeff(Wperp,J,M,0); b0=coeff(Wpar,J,M,0); a1=coeff(Wperp,J,M,1); b1=coeff(Wpar,J,M,1)
        w=max(w,abs(a0-1/(4*pi)),abs(b0-1/(4*pi)))
        pa1=mp.mpf(5)/(4*pi)*(J*(J+1)-3)*(J*(J+1)-3*M*M)/(J*(J+1)*(2*J-1)*(2*J+3))
        pb1=mp.mpf(5)/(4*pi)*(J*(J+1)-3*M*M)/((2*J-1)*(2*J+3))
        w=max(w,abs(a1-pa1),abs(b1-pb1))
ok&=report("7.3.159/160/161 a0=b0=1/4pi, a1, b1", w, tol=mp.mpf('1e-9'))

# ---- 7.3.168/169 special M ----
print("7.3.168/169 special-M forms")
def Pd(L,x): return float(mp.diff(lambda t: mp.legendre(L,t), x))
w=mp.mpf(0)
for J in [2,3,4]:
    # W^perp_{J0} = (2J+1)/(4pi J(J+1)) sin^2 [P_J']^2 ; W^par_{J0}=(2J+1)/4pi [P_J]^2
    w=max(w,abs(Wperp(J,0,TH)-float((2*J+1)/(4*pi*(J+1)*J)*Sn**2*Pd(J,C)**2)))
    w=max(w,abs(Wpar(J,0,TH)-float((2*J+1)/(4*pi)*mp.legendre(J,C)**2)))
    # W^par_{JJ} = (2J+1)!/(4pi 2^{2J}(J!)^2) sin^{2J}
    w=max(w,abs(Wpar(J,J,TH)-float(mp.factorial(2*J+1)/(4*pi*2**(2*J)*mp.factorial(J)**2)*Sn**(2*J))))
    # W^perp_{JJ} = (2J+1)!/(4pi 2^{2J}(J+1)!(J-1)!) sin^{2J-2}(1+cos^2)
    w=max(w,abs(Wperp(J,J,TH)-float(mp.factorial(2*J+1)/(4*pi*2**(2*J)*mp.factorial(J+1)*mp.factorial(J-1))*Sn**(2*J-2)*(1+C**2))))
ok&=report("7.3.168/169 W_{J0}, W_{JJ}", w, tol=mp.mpf('1e-10'))

# ---- Table 7.3.2 (corrected reading) ----
print("Table 7.3.2  W^perp / W^par (21 rows, OCR-corrected)")
def wp(J,M): return Wperp(J,M,TH)
def wl(J,M): return Wpar(J,M,TH)
perp={
 # (0,0): transverse VSH do not exist at J=0 (book lists 1/4pi); skip
 (1,0):3/(8*pi)*Sn**2, (1,1):3/(16*pi)*(1+C**2),
 (2,0):15/(8*pi)*Sn**2*C**2, (2,1):5/(16*pi)*(1-3*C**2+4*C**4), (2,2):5/(16*pi)*(1-C**4),
 (3,0):21/(64*pi)*Sn**2*(1-5*C**2)**2, (3,1):7/(256*pi)*(1+111*C**2-305*C**4+225*C**6),
 (3,2):35/(128*pi)*Sn**2*(1-2*C**2+9*C**4), (3,3):105/(256*pi)*Sn**4*(1+C**2),
 (4,0):45/(64*pi)*Sn**2*C**2*(3-7*C**2)**2, (4,1):9/(256*pi)*(9-153*C**2+855*C**4-1463*C**6+784*C**8),
 (4,2):9/(128*pi)*Sn**2*(1+50*C**2-175*C**4+196*C**6), (4,3):63/(256*pi)*Sn**4*(1+C**2+16*C**4),
 (4,4):63/(128*pi)*Sn**6*(1+C**2),
 (5,0):165/(512*pi)*Sn**2*(1-14*C**2+21*C**4)**2,
 (5,1):11/(1024*pi)*(1+813*C**2-7070*C**4+21378*C**6-26019*C**8+11025*C**10),
 (5,2):77/(256*pi)*Sn**2*(1-20*C**2+150*C**4-324*C**6+225*C**8),
 (5,3):231/(2048*pi)*Sn**4*(1+31*C**2-129*C**4+225*C**6),   # book "225 cos^8" -> cos^6 ?
 (5,4):231/(1024*pi)*Sn**6*(1+6*C**2+25*C**4), (5,5):1155/(2048*pi)*Sn**8*(1+C**2),
}
par={
 (0,0):1/(4*pi),
 (1,0):3/(4*pi)*C**2, (1,1):3/(8*pi)*Sn**2,
 (2,0):5/(16*pi)*(1-3*C**2)**2, (2,1):15/(8*pi)*Sn**2*C**2, (2,2):15/(32*pi)*Sn**4,
 (3,0):7/(16*pi)*C**2*(3-5*C**2)**2, (3,1):21/(64*pi)*Sn**2*(1-5*C**2)**2,
 (3,2):105/(32*pi)*Sn**4*C**2, (3,3):35/(64*pi)*Sn**6,
 (4,0):9/(256*pi)*(3-30*C**2+35*C**4)**2, (4,1):45/(64*pi)*Sn**2*C**2*(3-7*C**2)**2,
 (4,2):45/(128*pi)*Sn**4*(1-7*C**2)**2, (4,3):315/(64*pi)*Sn**6*C**2, (4,4):315/(512*pi)*Sn**8,
 (5,0):11/(256*pi)*C**2*(15-70*C**2+63*C**4)**2, (5,1):165/(512*pi)*Sn**2*(1-14*C**2+21*C**4)**2,
 (5,2):1155/(128*pi)*Sn**4*C**2*(1-3*C**2)**2, (5,3):385/(1024*pi)*Sn**6*(1-9*C**2)**2,  # book "365" is OCR of 385
 (5,4):3465/(512*pi)*Sn**8*C**2, (5,5):693/(1024*pi)*Sn**10,
}
wperp_bad=[]; wpar_bad=[]
for (J,M),expr in perp.items():
    if abs(wp(J,M)-expr)>1e-9: wperp_bad.append((J,M,float(wp(J,M)-expr)))
for (J,M),expr in par.items():
    if abs(wl(J,M)-expr)>1e-9: wpar_bad.append((J,M,float(wl(J,M)-expr)))
print("  W^perp mismatches:", wperp_bad if wperp_bad else "none")
print("  W^par  mismatches:", wpar_bad if wpar_bad else "none")
ok&= (not wperp_bad and not wpar_bad)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
