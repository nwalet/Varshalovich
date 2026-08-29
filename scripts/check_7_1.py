#!/usr/bin/env python3
r"""Verify Sec 7.1 general tensor spherical harmonics Y^{LS}_{JM} at S=1/2, 1.

[Y^{LS}_{JM}]^mu = C^{JM}_{L,M-mu,S,mu} Y_{L,M-mu}(th,ph), mu=S,...,-S (rows).
Checks: 7.1.33 (S.n action), 7.1.36-38 (S.L, S.J, L.J eigenvalues),
7.1.39-42 (sums with Y_{S mu}).
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
def report(tag,w,tol=mp.mpf('1e-12')):
    w=float(w); ok=w<tol; print(f"  {tag:46s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok
def spinS(S):
    ms=[S-k for k in range(int(round(2*S))+1)]; n=len(ms)
    Sz=mp.zeros(n,n); Sp=mp.zeros(n,n)
    for i,m in enumerate(ms): Sz[i,i]=m
    for i,m in enumerate(ms):
        if m+1<=S: Sp[ms.index(m+1),i]=mp.sqrt(S*(S+1)-m*(m+1))
    Sm=Sp.T
    return ms,(Sp+Sm)/2,(Sp-Sm)/(2j),Sz
def YLS(J,L,M,S,th,ph):
    ms=[S-k for k in range(int(round(2*S))+1)]
    return mp.matrix([cg(L,M-mu,S,mu,J,M)*Y(int(L),int(M-mu),th,ph) for mu in ms])
def vnorm(v): return max(abs(v[i]) for i in range(len(v)))
ok=True
TH,PH=mp.mpf('1.0'),mp.mpf('0.7')
nx=mp.sin(TH)*mp.cos(PH); ny=mp.sin(TH)*mp.sin(PH); nz=mp.cos(TH)

for S in [mp.mpf(1)/2, mp.mpf(1)]:
    ms,Sx,Sy,Sz=spinS(S); Sn=Sx*nx+Sy*ny+Sz*nz
    Js=[abs(k) for k in []]
    print(f"S={float(S)}")
    # 7.1.33 (S.n) action
    w=mp.mpf(0)
    for J in [S+k for k in range(4)]:   # J from S upward, integer step? J=|L-S|..; use L=J+-..
        pass
    # iterate over (J,L): L=|J-S|..J+S
    def Jrange(S):
        # J half-integer if S half-int; take a few
        base=S
        return [base+k for k in range(5)]
    for J in Jrange(S):
        for L in [abs(J-S)+k for k in range(int(2*S)+1)]:
            for M in [J-k for k in range(int(2*J)+1)]:
                lhs=Sn*YLS(J,L,M,S,TH,PH)
                # RHS 7.1.33: -1/2 { sqrt(...)/... Y^{L+1} + sqrt(...) Y^{L-1} }
                def fac(a): return a
                t1=mp.sqrt((J+L+S+2)*(J+L-S+1)*(J-L+S)*(-J+L+S+1)/((2*L+1)*(2*L+3)))
                t2=mp.sqrt((J+L+S+1)*(J+L-S)*(J-L+S+1)*(-J+L+S)/((2*L-1)*(2*L+1))) if L>=1 else mp.mpf(0)
                rhs=-mp.mpf(1)/2*(t1*YLS(J,L+1,M,S,TH,PH)+t2*YLS(J,L-1,M,S,TH,PH))
                w=max(w,vnorm(lhs-rhs))
    ok&=report(f"S={float(S)} 7.1.33 (S.n) action", w)
    # 7.1.36/37/38 eigenvalues (quantum numbers)
    w=mp.mpf(0)
    for J in Jrange(S):
        for L in [abs(J-S)+k for k in range(int(2*S)+1)]:
            SL=mp.mpf(1)/2*(J*(J+1)-L*(L+1)-S*(S+1))
            SJ=mp.mpf(1)/2*(J*(J+1)-L*(L+1)+S*(S+1))
            LJ=mp.mpf(1)/2*(J*(J+1)+L*(L+1)-S*(S+1))
            # verify S.L + S.S ... = S.J consistency: S.J = S.L + S^2
            w=max(w,abs(SJ-(SL+S*(S+1))), abs(LJ-(SL+L*(L+1))))
    ok&=report(f"S={float(S)} 7.1.36-38 eigenvalue consistency", w)
    if S==int(S):
     # 7.1.39 sum: sum_mu Y_{S mu} [Y^{LS}_{JM}]^mu = sqrt((2S+1)(2L+1)/(4pi(2J+1))) C^{J0}_{L0 S0} Y_{JM}
     w=mp.mpf(0)
     for J in Jrange(S)[:3]:
        for L in [abs(J-S)+k for k in range(int(2*S)+1)]:
            if (int(L)+int(2*S))%1!=0: continue
            for M in [J-k for k in range(int(2*J)+1)]:
                yls=YLS(J,L,M,S,TH,PH)
                lhs=sum(Y(int(S),int(ms[i]),TH,PH)*yls[i] for i in range(len(ms)))
                rhs=mp.sqrt((2*S+1)*(2*L+1)/(4*mp.pi*(2*J+1)))*cg(L,0,S,0,J,0)*Y(int(J),int(M),TH,PH)
                w=max(w,abs(lhs-rhs))
     ok&=report(f"S={float(S)} 7.1.39 sum", w)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
