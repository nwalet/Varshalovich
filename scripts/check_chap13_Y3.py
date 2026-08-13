from sympy.physics.wigner import gaunt
from sympy import sqrt, Rational, S, pi, N
from sympy.physics.wigner import clebsch_gordan as _CG0, wigner_6j as _6j0
def CG(*a):
    try: return _CG0(*a)
    except Exception: return S(0)
def wigner_6j(*a):
    try: return _6j0(*a)
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

# generic uncoupled apply of J_1beta to (ml,ms): returns list (ml',ms',amp) at same l
def Japply(l,ml,ms,beta,s):
    r=[]
    a=L_op(l,ml+beta,beta,l,ml)
    if a!=0: r.append((ml+beta,ms,a))
    b=S_op(s,ms+beta,beta,s,ms)
    if b!=0: r.append((ml,ms+beta,b))
    return r

# === 13.2.138  {n x J}_{k kappa} Y_{L nu} ===
def nJY_true(lp,s,Jp,Mp,k,kap,L,nu,l,J,M):
    tot=S(0)
    for al in(-1,0,1):
        be=kap-al
        if abs(be)>1: continue
        cg=CG(1,1,k,al,be,kap)
        if cg==0: continue
        for (mlk,msk,ck) in coupled(l,s,J,M):
            # Y_Lnu: orbital l-> l1 = lp? intermediate general; we sum over l1 by Y_op to arbitrary
            for l1 in range(abs(l-L),l+L+1):
                mY=mlk+nu
                aY=Y_op(l1,mY,L,nu,l,mlk)
                if aY==0: continue
                # J_1be on |l1 mY, s msk>
                for (mli,msi,aJ) in Japply(l1,mY,msk,be,s):
                    # n_1al: orbital l1->lp
                    an=n_op(lp,mli+al,al,l1,mli)
                    if an==0: continue
                    mlf=mli+al; msf=msi
                    for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                        if mlb==mlf and msb==msf: tot+=cb*ck*cg*aY*aJ*an
    return tot
def nJY_book(lp,s,Jp,Mp,k,kap,L,nu,l,J,M):
    if s!=Sp(s): pass
    pre=(-1)**(2*s+J-Jp+k+L)*sqrt((2*L+1)*(2*l+1)*(2*J+1)*(2*k+1)/(4*pi))
    tot=S(0)
    for L1 in range(abs(L-1),L+l+2):
        for J1 in [x*h for x in range(1,12)]:
            M1=M+nu
            if abs(M1)>J1: continue
            tot+=(sqrt(J1*(J1+1)*(2*L1+1))*(2*J1+1)
                *CG(J,L,J1,M,nu,M1)*CG(J1,k,Jp,M1,kap,Mp)
                *CG(L,l,L1,0,0,0)*CG(L1,1,lp,0,0,0)
                *wigner_6j(l,s,J,J1,L,L1)*wigner_6j(1,1,k,Jp,J1,J1)*wigner_6j(L1,s,J1,Jp,1,lp))
    return pre*tot
def Sp(x): return x
c138=[(1,h,3*h,h,1,0,1,0,1,3*h,h),(2,h,3*h,h,1,0,1,0,2,3*h,h),(1,h,h,h,2,0,1,0,1,3*h,h),(3,h,5*h,h,1,0,2,1,2,3*h,h)]
print("138 {n x J} Y (CG=C^{L1}_{L0,l0}):",all(close(nJY_true(*c),nJY_book(*c)) for c in c138))

# === 13.2.144  {L x J}_{k kappa} Y_{L nu}  (l'=l via L? no: L changes l) ===
def LJY_true(lp,s,Jp,Mp,k,kap,L,nu,l,J,M):
    tot=S(0)
    for al in(-1,0,1):
        be=kap-al
        if abs(be)>1: continue
        cg=CG(1,1,k,al,be,kap)
        if cg==0: continue
        for (mlk,msk,ck) in coupled(l,s,J,M):
            mY=mlk+nu
            aY=Y_op(lp,mY,L,nu,l,mlk)   # Y: l->l' (L diagonal after)
            if aY==0: continue
            for (mli,msi,aJ) in Japply(lp,mY,msk,be,s):
                aL=L_op(lp,mli+al,al,lp,mli)
                if aL==0: continue
                mlf=mli+al; msf=msi
                for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                    if mlb==mlf and msb==msf: tot+=cb*ck*cg*aY*aJ*aL
    return tot
def LJY_book(lp,s,Jp,Mp,k,kap,L,nu,l,J,M):
    pre=(-1)**(l+lp+J+Jp+k+1)*sqrt((2*L+1)*(2*l+1)*(2*J+1)*(2*k+1)*lp*(lp+1)*(2*lp+1)/(4*pi))*CG(L,l,lp,0,0,0)
    tot=S(0)
    for J1 in [x*h for x in range(1,12)]:
        M1=M+nu
        if abs(M1)>J1: continue
        tot+=((2*J1+1)*sqrt(J1*(J1+1))*CG(J,L,J1,M,nu,M1)*CG(J1,k,Jp,M1,kap,Mp)
            *wigner_6j(l,s,J,J1,L,lp)*wigner_6j(1,1,k,Jp,J1,J1)*wigner_6j(lp,s,J1,Jp,1,lp))
    return pre*tot
c144=[(1,h,3*h,h,1,0,1,0,2,3*h,h),(2,h,3*h,h,1,0,1,0,2,3*h,h),(3,h,5*h,h,1,0,2,1,2,3*h,h),(1,h,h,h,2,0,1,0,1,3*h,h)]
print("144 {L x J} Y (CG=C^{l'0}_{L0,l0}):",all(close(LJY_true(*c),LJY_book(*c)) for c in c144))
