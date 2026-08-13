"""Verify Y_L coupled tensor/commutator MEs, Chap13 sec 13.2.8 (123,125,126,127)."""
from sympy.physics.wigner import clebsch_gordan as CG, gaunt, wigner_6j, wigner_9j
from sympy import sqrt, Rational, simplify, S, pi, N
def Y_op(lp,mp,L,nu,l,m):
    if lp<0 or L<0: return S(0)
    return (-1)**mp*gaunt(lp,L,l,-mp,nu,m)
def L_op(lp,mlp,mu,l,ml):
    if lp!=l or mlp!=ml+mu: return S(0)
    return sqrt(l*(l+1))*CG(l,1,l,ml,mu,mlp)
def S_op(sp,msp,mu,s,ms):
    if sp!=s or msp!=ms+mu: return S(0)
    return sqrt(s*(s+1))*CG(s,1,s,ms,mu,msp)
h=Rational(1,2)
def coupled(l,s,J,M):
    out=[]
    for i in range(2*l+1):
        ml=-l+i; ms=M-ml
        if abs(ms)<=s:
            c=CG(l,s,J,ml,ms,M)
            if c!=0: out.append((ml,ms,c))
    return out
def J_action(l,ml,ms,nu,s):
    r=[]
    a=L_op(l,ml+nu,nu,l,ml)
    if a!=0: r.append((ml+nu,ms,a))
    b=S_op(s,ms+nu,nu,s,ms)
    if b!=0: r.append((ml,ms+nu,b))
    return r
def close(a,b): return abs(complex(N(a-b)))<1e-9

# ---- 13.2.123 {Y_L x S_1}_{L'nu'} ----
def YS_true(lp,s,Jp,Mp,Lp,vp,L,l,J,M):
    tot=S(0)
    for mu in(-1,0,1):
        nu=vp-mu
        cg=CG(L,1,Lp,nu,mu,vp)
        if cg==0: continue
        for (mlk,msk,ck) in coupled(l,s,J,M):
            for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                amp=Y_op(lp,mlb,L,nu,l,mlk)*S_op(s,msb,mu,s,msk)
                if amp!=0: tot+=cb*ck*cg*amp
    return tot
def YS_book(lp,s,Jp,Mp,Lp,vp,L,l,J,M):
    return (sqrt((2*J+1)*(2*Lp+1)*(2*L+1)*(2*l+1)*s*(s+1)*(2*s+1)/(4*pi))
        *CG(l,L,lp,0,0,0)*CG(J,Lp,Jp,M,vp,Mp)*wigner_9j(l,L,lp,s,1,s,J,Lp,Jp))
c123=[(2,h,5*h,h,2,0,2,2,3*h,h),(1,h,3*h,h,1,0,1,1,h,h),(2,h,3*h,3*h,2,1,1,2,h,h),(0,h,h,h,1,0,1,1,h,h)]
print("123 {Y x S}:",all(close(YS_true(*c),YS_book(*c)) for c in c123))

# ---- 13.2.125 {J_1 x Y_L}_{L'nu'} ----
def JYt_true(lp,s,Jp,Mp,Lp,vp,L,l,J,M):
    tot=S(0)
    for mu in(-1,0,1):
        nu=vp-mu
        cg=CG(1,L,Lp,mu,nu,vp)
        if cg==0: continue
        for (mlk,msk,ck) in coupled(l,s,J,M):
            mli=mlk+nu
            aY=Y_op(lp,mli,L,nu,l,mlk)
            if aY==0: continue
            for (mlf,msf,ampJ) in J_action(lp,mli,msk,mu,s):
                for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                    if mlb==mlf and msb==msf: tot+=cb*ck*cg*aY*ampJ
    return tot
def JYt_book(lp,s,Jp,Mp,Lp,vp,L,l,J,M):   # confirmed: phase (-1)^{l+s+L'-J'}, CG C^{l'0}_{l0,L0}
    return ((-1)**(l+s+Lp-Jp)*sqrt((2*L+1)*(2*Lp+1)*(2*l+1)*(2*J+1)*Jp*(Jp+1)*(2*Jp+1)/(4*pi))
        *CG(J,Lp,Jp,M,vp,Mp)*wigner_6j(L,1,Lp,Jp,J,Jp)*CG(l,L,lp,0,0,0)*wigner_6j(l,s,J,Jp,L,lp))
c125=[(1,h,3*h,h,1,0,1,2,3*h,h),(3,h,5*h,h,1,0,1,2,3*h,h),(1,h,h,h,2,0,2,2,3*h,h),(1,h,3*h,3*h,2,1,2,3,5*h,h)]
print("125 {J x Y}:",all(close(JYt_true(*c),JYt_book(*c)) for c in c125))

# ---- 13.2.126 <..|[J_1mu,Y_Lnu]|..> = [L,Y] ----
def commJY_true(lp,s,Jp,Mp,mu,L,nu,l,J,M):
    tot=S(0)
    for (mlk,msk,ck) in coupled(l,s,J,M):
        # L_1mu Y_Lnu
        mli=mlk+nu
        a=Y_op(lp,mli,L,nu,l,mlk)*L_op(lp,mli+mu,mu,lp,mli)
        # Y_Lnu L_1mu (L diagonal on ket l): L then Y
        b=L_op(l,mlk+mu,mu,l,mlk)*Y_op(lp,mlk+mu+nu,L,nu,l,mlk+mu)
        for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
            if msb==msk:
                if mlb==mli+mu: tot+=cb*ck*a
                if mlb==mlk+mu+nu: tot-=cb*ck*b
    return tot
def commJY_book(lp,s,Jp,Mp,mu,L,nu,l,J,M):   # confirmed: phase (-1)^{J+l'+s+L}, CG C^{l'0}_{l0,L0}
    return ((-1)**(J+lp+s+L)*sqrt(L*(L+1)*(2*L+1)*(2*l+1)*(2*J+1)/(4*pi))
        *CG(l,L,lp,0,0,0)*wigner_6j(l,s,J,Jp,L,lp)*CG(L,1,L,nu,mu,mu+nu)*CG(J,L,Jp,M,mu+nu,Mp))
c126=[(1,h,3*h,h,0,1,0,2,3*h,h),(3,h,5*h,h,1,1,-1,2,3*h,h),(1,h,h,h,0,2,0,2,3*h,h),(3,h,7*h,h,0,2,0,1,h,h)]
print("126 [J,Y] direct:",all(close(commJY_true(*c),commJY_book(*c)) for c in c126))

# ---- 13.2.127 irreducible commutator R^{1L}_{L'nu'}(J,Y): delta_{LL'}, comp nu'=nu ----
def commJYirr_true(lp,s,Jp,Mp,Lp,vp,L,l,J,M):
    # R^{1L}_{L' nu'} = -sqrt(L(L+1)) Y_{L,nu'} delta_{L'L}; ME via [L,Y] projected
    if Lp!=L: return S(0)
    tot=S(0)
    for (mlk,msk,ck) in coupled(l,s,J,M):
        for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
            if msb==msk:
                amp=-sqrt(L*(L+1))*Y_op(lp,mlb,L,vp,l,mlk)
                if amp!=0: tot+=cb*ck*amp
    return tot
def commJYirr_book(lp,s,Jp,Mp,Lp,vp,L,l,J,M):   # confirmed: phase (-1)^{J+l'+s+L+1}, CG C^{l'0}_{l0,L0}
    if Lp!=L: return S(0)
    return ((-1)**(J+lp+s+L+1)*sqrt(L*(L+1)*(2*L+1)*(2*l+1)*(2*J+1)/(4*pi))
        *CG(l,L,lp,0,0,0)*wigner_6j(l,s,J,Jp,L,lp)*CG(J,L,Jp,M,vp,Mp))
c127=[(1,h,3*h,h,1,0,1,2,3*h,h),(3,h,5*h,h,1,-1,1,2,3*h,h),(1,h,h,h,2,0,2,2,3*h,h),(3,h,7*h,h,2,0,2,1,h,h)]
print("127 [J,Y] irreducible:",all(close(commJYirr_true(*c),commJYirr_book(*c)) for c in c127))
