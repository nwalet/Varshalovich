"""Verify spin angular momentum S matrix elements, Chap13 sec 13.2.7."""
from sympy.physics.wigner import clebsch_gordan as CG, gaunt, wigner_6j, wigner_9j
from sympy import sqrt, Rational, simplify, S, pi

def L_op(lp,mlp,mu,l,ml):
    if lp!=l or mlp!=ml+mu: return S(0)
    return sqrt(l*(l+1))*CG(l,1,l,ml,mu,mlp)
def S_op(sp,msp,mu,s,ms):
    if sp!=s or msp!=ms+mu: return S(0)
    return sqrt(s*(s+1))*CG(s,1,s,ms,mu,msp)
def n_op(lp,mlp,mu,l,ml):
    if lp<0: return S(0)
    return sqrt(4*pi/3)*(-1)**mlp*gaunt(lp,1,l,-mlp,mu,ml)

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

# ---- 13.2.96 reduced ME <l s J'||S||l s J> ----
print("=== 13.2.96 reduced ME of S in (lsJ) ===")
def S_me_coupled(l,s,Jp,Mp,mu,J,M):
    tot=S(0)
    for (mlk,msk,ck) in coupled(l,s,J,M):
        for (mlb,msb,cb) in coupled(l,s,Jp,Mp):
            if mlb==mlk:
                tot+=cb*ck*S_op(s,msb,mu,s,msk)
    return simplify(tot)
def redS_book(l,s,Jp,J):
    return simplify((-1)**(l+s+Jp+1)*sqrt(s*(s+1)*(2*s+1)*(2*J+1)*(2*Jp+1))*wigner_6j(s,l,J,Jp,1,s))
def redS_from_me(l,s,Jp,Mp,mu,J,M):
    cg=CG(J,1,Jp,M,mu,Mp)
    if cg==0: return None
    return simplify(S_me_coupled(l,s,Jp,Mp,mu,J,M)*sqrt(2*Jp+1)/cg)
ok=True
for (l,J,Jp) in [(1,h,h),(1,3*h,3*h),(1,3*h,h),(2,3*h,5*h),(2,5*h,3*h),(0,h,h)]:
    for mu in (0,1,-1):
        M=h; Mp=M+mu
        if abs(Mp)>Jp or abs(M)>J: continue
        r=redS_from_me(l,h,Jp,Mp,mu,J,M)
        if r is None: continue
        if simplify(r-redS_book(l,h,Jp,J))!=0: ok=False
print("  ",ok)

# ---- generic orbital(x)spin tensor product ----
def tps(Afun,Bfun,lp,s,Jp,Mp,k,kap,l,J,M):
    tot=S(0)
    for mu in(-1,0,1):
        nu=kap-mu
        if abs(nu)>1: continue
        cg=CG(1,1,k,mu,nu,kap)
        if cg==0: continue
        for (mlk,msk,ck) in coupled(l,s,J,M):
            for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                amp=Afun(lp,mlb,mu,l,mlk)*Bfun(s,msb,nu,s,msk)
                if amp!=0: tot+=cb*ck*cg*amp
    return simplify(tot)

print("=== 13.2.97 {n x S}_{k kappa} ===")
def nS_book(lp,s,Jp,Mp,k,kap,l,J,M):
    return simplify(sqrt((2*k+1)*(2*J+1)*(2*l+1)*(2*s+1)*s*(s+1))
        *CG(l,1,lp,0,0,0)*CG(J,k,Jp,M,kap,Mp)*wigner_9j(l,1,lp,s,1,s,J,k,Jp))
print("=== 13.2.99 {L x S}_{k kappa} ===")
def LS_book(lp,s,Jp,Mp,k,kap,l,J,M):
    return simplify(sqrt((2*k+1)*(2*J+1)*(2*l+1)*(2*s+1)*l*(l+1)*s*(s+1))
        *CG(J,k,Jp,M,kap,Mp)*wigner_9j(l,1,l,s,1,s,J,k,Jp))
cc=[(1,h,h,h,1,0,1,h,h),(1,h,3*h,h,1,0,1,3*h,h),(2,h,5*h,h,1,0,2,3*h,h),
    (1,h,3*h,3*h,2,1,1,3*h,h),(2,h,3*h,h,2,0,2,3*h,h),(0,h,h,h,1,0,0,h,h)]
print("  97:",all(simplify(tps(n_op,S_op,*c)-nS_book(*c))==0 for c in cc))
cc2=[c for c in cc if c[0]==c[6]]   # L x S needs l'=l
print("  99:",all(simplify(tps(L_op,S_op,*c)-LS_book(*c))==0 for c in cc2))

# ---- {S x J} and {J x S} ----
print("=== 13.2.100 {S x J} and 13.2.101 {J x S} ===")
def SJ_true(lp,s,Jp,Mp,k,kap,l,J,M):
    tot=S(0)
    if lp!=l: return tot
    for mu in(-1,0,1):
        nu=kap-mu
        if abs(nu)>1: continue
        cg=CG(1,1,k,mu,nu,kap)
        if cg==0: continue
        for (mlk,msk,ck) in coupled(l,s,J,M):
            for (mli,msi,ampJ) in J_action(l,mlk,msk,nu,s):
                aS=S_op(s,msi+mu,mu,s,msi)
                if aS==0: continue
                mlf=mli; msf=msi+mu
                for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                    if mlb==mlf and msb==msf: tot+=cb*ck*cg*ampJ*aS
    return simplify(tot)
def JS_true(lp,s,Jp,Mp,k,kap,l,J,M):
    tot=S(0)
    if lp!=l: return tot
    for mu in(-1,0,1):
        nu=kap-mu
        if abs(nu)>1: continue
        cg=CG(1,1,k,mu,nu,kap)
        if cg==0: continue
        for (mlk,msk,ck) in coupled(l,s,J,M):
            aS=S_op(s,msk+nu,nu,s,msk)
            if aS==0: continue
            mli=mlk; msi=msk+nu
            for (mlf,msf,ampJ) in J_action(l,mli,msi,mu,s):
                for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                    if mlb==mlf and msb==msf: tot+=cb*ck*cg*aS*ampJ
    return simplify(tot)
def SJ_book(lp,s,Jp,Mp,k,kap,l,J,M):
    return simplify((-1)**(J+k-l-s+1)*(2*J+1)*sqrt((2*k+1)*s*(s+1)*(2*s+1)*J*(J+1))
        *wigner_6j(1,1,k,Jp,J,J)*wigner_6j(s,l,J,Jp,1,s)*CG(J,k,Jp,M,kap,Mp))
def JS_book(lp,s,Jp,Mp,k,kap,l,J,M):
    return simplify((-1)**(J+k-l-s+1)*sqrt((2*k+1)*(2*J+1)*Jp*(Jp+1)*(2*Jp+1))
        *sqrt(s*(s+1)*(2*s+1))*wigner_6j(1,1,k,Jp,J,Jp)*wigner_6j(s,l,J,Jp,1,s)*CG(J,k,Jp,M,kap,Mp))
print("  100:",all(simplify(SJ_true(*c)-SJ_book(*c))==0 for c in cc2))
print("  101:",all(simplify(JS_true(*c)-JS_book(*c))==0 for c in cc2))
