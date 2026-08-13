from sympy.physics.wigner import gaunt
from sympy import sqrt, Rational, S, pi, N
from sympy.physics.wigner import clebsch_gordan as _CG0, wigner_6j as _6j0, wigner_9j as _9j0
def CG(*a):
    try: return _CG0(*a)
    except Exception: return S(0)
def w6(*a):
    try: return _6j0(*a)
    except Exception: return S(0)
def w9(*a):
    try: return _9j0(*a)
    except Exception: return S(0)
h=Rational(1,2)
def Y_op(lp,mp,L,nu,l,m):
    if lp<0 or L<0: return S(0)
    return (-1)**mp*gaunt(lp,L,l,-mp,nu,m)
def n_op(lp,mlp,mu,l,ml):
    if lp<0: return S(0)
    return sqrt(4*pi/3)*(-1)**mlp*gaunt(lp,1,l,-mlp,mu,ml)
def L_op(lp,mlp,mu,l,ml):
    if lp!=l or mlp!=ml+mu: return S(0)
    return sqrt(l*(l+1))*CG(l,1,l,ml,mu,mlp)
def S_op(sp,msp,mu,s,ms):
    if sp!=s or msp!=ms+mu: return S(0)
    return sqrt(s*(s+1))*CG(s,1,s,ms,mu,msp)
def coupled(l,s,J,M):
    out=[]
    for i in range(2*l+1):
        ml=-l+i; ms=M-ml
        if abs(ms)<=s:
            c=CG(l,s,J,ml,ms,M)
            if c!=0: out.append((ml,ms,c))
    return out
def close(a,b): return abs(complex(N(a-b)))<1e-9
Js=[x*h for x in range(0,14)]

# ===== 13.2.130  {n x L}_{k kappa} Y_{L nu}  (orbital) =====
def nLY_true(lp,mp,k,kap,L,nu,l,m):
    tot=S(0)
    for al in(-1,0,1):
        be=kap-al
        if abs(be)>1: continue
        cg=CG(1,1,k,al,be,kap)
        if cg==0: continue
        for lam in range(abs(l-L),l+L+1):
            mY=m+nu
            aY=Y_op(lam,mY,L,nu,l,m)
            if aY==0: continue
            aL=L_op(lam,mY+be,be,lam,mY)   # L diagonal
            if aL==0: continue
            an=n_op(lp,mY+be+al,al,lam,mY+be)
            if an==0: continue
            tot+=cg*aY*aL*an
    return tot
def nLY_book(lp,mp,k,kap,L,nu,l,m):
    tot=S(0)
    for Lp in range(abs(l-L),l+L+1):
        vp=m+nu; mpp=mp
        tot+=((-1)**(lp+Lp+k)*sqrt((2*L+1)*(2*l+1)*(2*k+1)*Lp*(Lp+1)*(2*Lp+1)/(4*pi*(2*lp+1)))
            *CG(l,L,Lp,m,nu,vp)*CG(Lp,k,lp,vp,kap,mp)*CG(l,L,Lp,0,0,0)*CG(Lp,1,lp,0,0,0)
            *w6(1,1,k,lp,Lp,Lp))
    return tot
c130=[(2,1,1,0,1,0,1,1),(1,0,2,1,1,0,2,0),(3,1,1,0,2,1,2,0),(2,0,1,-1,1,0,2,1)]
print("130 {n x L} Y:",all(close(nLY_true(*c),nLY_book(*c)) for c in c130))

# ===== 13.2.134  {n x S}_{k kappa} Y_{L nu}  (coupled) =====
def nSY_true(lp,s,Jp,Mp,k,kap,L,nu,l,J,M):
    tot=S(0)
    for al in(-1,0,1):
        be=kap-al
        if abs(be)>1: continue
        cg=CG(1,1,k,al,be,kap)
        if cg==0: continue
        for (mlk,msk,ck) in coupled(l,s,J,M):
            for lam in range(abs(l-L),l+L+1):
                mY=mlk+nu
                aY=Y_op(lam,mY,L,nu,l,mlk)
                if aY==0: continue
                aS=S_op(s,msk+be,be,s,msk)
                if aS==0: continue
                an=n_op(lp,mY+al,al,lam,mY)
                if an==0: continue
                mlf=mY+al; msf=msk+be
                for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                    if mlb==mlf and msb==msf: tot+=cb*ck*cg*aY*aS*an
    return tot
def nSY_book(lp,s,Jp,Mp,k,kap,L,nu,l,J,M,cg3):
    pre=(-1)**(l+s+J)*sqrt((2*L+1)*(2*l+1)*(2*J+1)*(2*k+1)*s*(s+1)*(2*s+1)/(4*pi))
    tot=S(0)
    for L1 in range(abs(l-L),l+L+1):
        for J1 in Js:
            M1=M+nu
            if abs(M1)>J1: continue
            tot+=(sqrt((2*J1+1)*(2*L1+1))*CG(J,L,J1,M,nu,M1)*CG(J1,k,Jp,M1,kap,Mp)
                *cg3(L,l,L1)*CG(L1,1,lp,0,0,0)
                *w6(l,s,J,J1,L,L1)*w9(L1,1,lp,s,1,s,J1,k,Jp))
    return pre*tot
c134=[(1,h,3*h,h,1,0,1,0,1,3*h,h),(2,h,3*h,h,1,0,1,0,2,3*h,h),(3,h,5*h,h,1,0,2,1,2,3*h,h),(1,h,h,h,2,0,1,0,1,3*h,h)]
print("134 (CG3=C^{L1}_{L0,l0}):",all(close(nSY_true(*c),nSY_book(*c,lambda L,l,L1:CG(L,l,L1,0,0,0))) for c in c134))

# ===== 13.2.142  {L x S}_{k kappa} Y_{L nu}  (coupled) =====
def LSY_true(lp,s,Jp,Mp,k,kap,L,nu,l,J,M):
    tot=S(0)
    for al in(-1,0,1):
        be=kap-al
        if abs(be)>1: continue
        cg=CG(1,1,k,al,be,kap)
        if cg==0: continue
        for (mlk,msk,ck) in coupled(l,s,J,M):
            mY=mlk+nu
            aY=Y_op(lp,mY,L,nu,l,mlk)   # Y: l->l', then L diagonal on l'
            if aY==0: continue
            aL=L_op(lp,mY+al,al,lp,mY)
            if aL==0: continue
            aS=S_op(s,msk+be,be,s,msk)
            if aS==0: continue
            mlf=mY+al; msf=msk+be
            for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                if mlb==mlf and msb==msf: tot+=cb*ck*cg*aY*aL*aS
    return tot
def LSY_book(lp,s,Jp,Mp,k,kap,L,nu,l,J,M,cg3,withpre):
    pre=(-1)**(l+s+Jp)
    if withpre: pre*=sqrt((2*L+1)*(2*l+1)*(2*J+1)*(2*k+1)*s*(s+1)*(2*s+1)*lp*(lp+1)*(2*lp+1)/(4*pi))
    tot=S(0)
    for J1 in Js:
        M1=M+nu
        if abs(M1)>J1: continue
        tot+=(sqrt(2*J1+1)*CG(J,L,J1,M,nu,M1)*CG(J1,k,Jp,M1,kap,Mp)
            *w6(l,s,J,J1,L,lp)*w9(lp,1,lp,s,1,s,J1,k,Jp))
    return pre*cg3*tot
c142=[(1,h,3*h,h,1,0,1,0,2,3*h,h),(2,h,3*h,h,1,0,1,0,2,3*h,h),(3,h,5*h,h,1,0,2,1,2,3*h,h),(1,h,h,h,2,0,1,0,1,3*h,h)]
print("142 (CG3=C^{l'0}_{L0,l0}, WITH sqrt prefactor):",
      all(close(LSY_true(*c),LSY_book(*c,CG(c[6],c[8],c[0],0,0,0),True)) for c in c142))
print("142 (CG3=C^{l'0}_{L0,10},  no prefactor):",
      all(close(LSY_true(*c),LSY_book(*c,CG(c[6],1,c[0],0,0,0),False)) for c in c142))

# ===== 13.2.146  {S x J}_{k kappa} Y_{L nu}  (coupled) =====
def SJY_true(lp,s,Jp,Mp,k,kap,L,nu,l,J,M):
    tot=S(0)
    def Japply(l1,ml,ms,be):
        r=[]
        a=L_op(l1,ml+be,be,l1,ml)
        if a!=0: r.append((ml+be,ms,a))
        b=S_op(s,ms+be,be,s,ms)
        if b!=0: r.append((ml,ms+be,b))
        return r
    for al in(-1,0,1):
        be=kap-al
        if abs(be)>1: continue
        cg=CG(1,1,k,al,be,kap)
        if cg==0: continue
        for (mlk,msk,ck) in coupled(l,s,J,M):
            mY=mlk+nu
            aY=Y_op(lp,mY,L,nu,l,mlk)   # Y: l->l'
            if aY==0: continue
            for (mli,msi,aJ) in Japply(lp,mY,msk,be):  # J_1be
                aS=S_op(s,msi+al,al,s,msi)             # S_1al (spin)
                if aS==0: continue
                mlf=mli; msf=msi+al
                for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                    if mlb==mlf and msb==msf: tot+=cb*ck*cg*aY*aJ*aS
    return tot
def SJY_book(lp,s,Jp,Mp,k,kap,L,nu,l,J,M):
    pre=(-1)**(J+L+k+1)*sqrt((2*L+1)*(2*l+1)*(2*J+1)*(2*k+1)*s*(s+1)*(2*s+1)/(4*pi))*CG(L,l,lp,0,0,0)
    tot=S(0)
    for J1 in Js:
        M1=M+nu
        if abs(M1)>J1: continue
        tot+=(CG(J,L,J1,M,nu,M1)*CG(J1,k,Jp,M1,kap,Mp)*sqrt(J1*(J1+1))*(2*J1+1)*(-1)**J1
            *w6(l,s,J,J1,L,lp)*w6(1,1,k,Jp,J1,J1)*w6(s,lp,J1,Jp,1,s))
    return pre*tot
c146=[(1,h,3*h,h,1,0,1,0,2,3*h,h),(2,h,3*h,h,1,0,1,0,2,3*h,h),(3,h,5*h,h,1,0,2,1,2,3*h,h),(1,h,h,h,2,0,1,0,1,3*h,h)]
print("146 {S x J} Y (CG=C^{l'0}_{L0,l0}):",all(close(SJY_true(*c),SJY_book(*c)) for c in c146))

# ===== F-coupled siblings: reuse direct-Y true, couple with C^{F phi}_{k kappa, L nu} =====
def Fcouple(true_fn,lp,s,Jp,Mp,k,F,phi,L,l,J,M):
    tot=S(0)
    for kap in range(-k,k+1):
        nu=phi-kap
        if abs(nu)>L: continue
        c=CG(k,L,F,kap,nu,phi)
        if c==0: continue
        tot+=c*true_fn(lp,s,Jp,Mp,k,kap,L,nu,l,J,M)
    return tot
# 13.2.135 {{n x S}_k x Y}_F
def nSY135_book(lp,s,Jp,Mp,k,F,phi,L,l,J,M):
    pre=(-1)**(l+s-Jp-F)*CG(J,F,Jp,M,phi,Mp)*sqrt((2*L+1)*(2*J+1)*(2*l+1)*(2*k+1)*(2*F+1)*s*(s+1)*(2*s+1)/(4*pi))
    tot=S(0)
    for L1 in range(abs(l-L),l+L+1):
        for J1 in Js:
            tot+=(sqrt(2*L1+1)*(2*J1+1)*CG(L,l,L1,0,0,0)*CG(L1,1,lp,0,0,0)
                *w6(l,s,J,J1,L,L1)*w6(L,k,F,Jp,J,J1)*w9(L1,1,lp,s,1,s,J1,k,Jp))
    return pre*tot
c135=[(1,h,3*h,h,1,1,0,1,1,3*h,h),(2,h,3*h,h,1,1,0,1,2,3*h,h),(3,h,5*h,h,1,2,1,2,2,5*h,h)]
print("135 {{n x S}_k x Y}_F:",all(close(Fcouple(nSY_true,*c),nSY135_book(*c)) for c in c135))

# 13.2.143 {{L x S}_k x Y}_F
def LSY143_book(lp,s,Jp,Mp,k,F,phi,L,l,J,M):
    pre=((-1)**(l+s-Jp-F)*CG(J,F,Jp,M,phi,Mp)*CG(L,l,lp,0,0,0)
        *sqrt((2*L+1)*(2*l+1)*(2*k+1)*(2*F+1)*s*(s+1)*(2*s+1)*lp*(lp+1)*(2*lp+1)*(2*J+1)/(4*pi)))
    tot=S(0)
    for J1 in Js:
        tot+=(2*J1+1)*w6(l,s,J,J1,L,lp)*w6(L,k,F,Jp,J,J1)*w9(lp,1,lp,s,1,s,J1,k,Jp)
    return pre*tot
c143=[(1,h,3*h,h,1,1,0,1,2,3*h,h),(2,h,3*h,h,1,1,0,1,2,3*h,h),(3,h,5*h,h,1,2,1,2,2,5*h,h)]
print("143 {{L x S}_k x Y}_F:",all(close(Fcouple(LSY_true,*c),LSY143_book(*c)) for c in c143))
