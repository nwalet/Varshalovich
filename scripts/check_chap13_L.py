"""Verify orbital angular momentum L matrix elements, Chap13 sec 13.2.6."""
from sympy.physics.wigner import clebsch_gordan as CG, gaunt, wigner_6j
from sympy import sqrt, Rational, simplify, S, pi

# ---------- single-space spherical operators (uncoupled basis) ----------
def L_op(lp,mlp,mu,l,ml):        # <l ml'|L_1mu|l ml>, diagonal in l
    if lp!=l or mlp!=ml+mu: return S(0)
    return sqrt(l*(l+1))*CG(l,1,l,ml,mu,mlp)
def S_op(sp,msp,mu,s,ms):
    if sp!=s or msp!=ms+mu: return S(0)
    return sqrt(s*(s+1))*CG(s,1,s,ms,mu,msp)
def n_op(lp,mlp,mu,l,ml):        # <l' ml'|n_1mu|l ml>  via Gaunt
    if lp<0: return S(0)
    return sqrt(4*pi/3)*(-1)**mlp*gaunt(lp,1,l,-mlp,mu,ml)

# =========================================================================
print("=== 13.2.73  {n x L}_{k kappa}  (orbital only) ===")
def nL_true(lp,mp,k,kap,l,m):
    t=S(0)
    for mu in(-1,0,1):
        nu=kap-mu
        if abs(nu)>1: continue
        t+=CG(1,1,k,mu,nu,kap)*n_op(lp,mp,mu,l,m+nu)*L_op(l,m+nu,nu,l,m)
    return simplify(t)
def nL_book(lp,mp,k,kap,l,m):
    return simplify((-1)**(lp+l+k)*CG(l,k,lp,m,kap,mp)
        *sqrt((2*k+1)*l*(l+1)/S(2*lp+1))*(2*l+1)*CG(l,1,lp,0,0,0)*wigner_6j(1,1,k,lp,l,l))
cases=[(2,1,1,0,2,0),(2,1,2,1,1,0),(3,2,1,1,2,1),(1,0,2,-1,2,1),(2,0,1,-1,2,1),(0,0,1,0,1,0)]
print("  ",all(simplify(nL_true(*c)-nL_book(*c))==0 for c in cases))

print("=== 13.2.74  {L x n}_{k kappa}  (orbital only) ===")
def Ln_true(lp,mp,k,kap,l,m):
    t=S(0)
    for mu in(-1,0,1):
        nu=kap-mu
        if abs(nu)>1: continue
        t+=CG(1,1,k,mu,nu,kap)*L_op(lp,mp,mu,lp,m+nu)*n_op(lp,m+nu,nu,l,m)
    return simplify(t)
def Ln_book(lp,mp,k,kap,l,m):
    return simplify((-1)**(lp+l+k)*CG(l,k,lp,m,kap,mp)
        *sqrt((2*k+1)*lp*(lp+1)*(2*l+1))*CG(l,1,lp,0,0,0)*wigner_6j(1,1,k,lp,l,lp))
print("  ",all(simplify(Ln_true(*c)-Ln_book(*c))==0 for c in cases))

# =========================================================================
# coupled (l s J) basis, s=1/2
h=Rational(1,2)
def coupled(l,s,J,M):
    out=[]
    for i in range(2*l+1):
        ml=-l+i
        ms=M-ml
        if abs(ms)<=s:
            c=CG(l,s,J,ml,ms,M)
            if c!=0: out.append((ml,ms,c))
    return out
def J_action(l,ml,ms,nu,s):      # J_1nu on |l ml>|s ms> -> list (mli,msi,amp)
    r=[]
    a=L_op(l,ml+nu,nu,l,ml)
    if a!=0: r.append((ml+nu,ms,a))
    b=S_op(s,ms+nu,nu,s,ms)
    if b!=0: r.append((ml,ms+nu,b))
    return r

print("=== 13.2.68  reduced ME <l s J'||L||l s J> ===")
def L_me_coupled(l,s,Jp,Mp,mu,J,M):   # <l s J' M'|L_1mu|l s J M>
    tot=S(0)
    for (mlk,msk,ck) in coupled(l,s,J,M):
        for (mlb,msb,cb) in coupled(l,s,Jp,Mp):
            if msb==msk:
                tot+=cb*ck*L_op(l,mlb,mu,l,mlk)
    return simplify(tot)
def red_book(l,s,Jp,J):
    return simplify((-1)**(s+J+l+1)*sqrt(2*J+1)*sqrt(2*Jp+1)*sqrt(2*l+1)*sqrt(l*(l+1))*wigner_6j(l,s,J,Jp,1,l))
def red_from_me(l,s,Jp,Mp,mu,J,M):    # invert Wigner-Eckart
    cg=CG(J,1,Jp,M,mu,Mp)
    if cg==0: return None
    return simplify(L_me_coupled(l,s,Jp,Mp,mu,J,M)*sqrt(2*Jp+1)/cg)
okred=True
for (l,J,Jp) in [(1,h,h),(1,3*h,3*h),(1,3*h,h),(2,3*h,5*h),(2,5*h,5*h),(1,h,3*h)]:
    for M in [Rational(1,2)]:
        for mu in (0,1,-1):
            Mp=M+mu
            if abs(Mp)>Jp: continue
            r=red_from_me(l,h,Jp,Mp,mu,J,M)
            if r is None: continue
            if simplify(r-red_book(l,h,Jp,J))!=0: okred=False
print("  ",okred)

print("=== 13.2.86  {L x J}_{k kappa}  and  13.2.87  {J x L}_{k kappa} ===")
def LJ_true(lp,s,Jp,Mp,k,kap,l,J,M):
    tot=S(0)
    for mu in(-1,0,1):
        for nu in(-1,0,1):
            if mu+nu!=kap: continue
            cg=CG(1,1,k,mu,nu,kap)
            if cg==0: continue
            for (mlk,msk,ck) in coupled(l,s,J,M):
                for (mli,msi,ampJ) in J_action(l,mlk,msk,nu,s):
                    if lp!=l: continue
                    mlf=mli+mu
                    for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                        if mlb==mlf and msb==msi:
                            tot+=cb*ck*cg*ampJ*L_op(l,mlf,mu,l,mli)
    return simplify(tot)
def JL_true(lp,s,Jp,Mp,k,kap,l,J,M):
    tot=S(0)
    for mu in(-1,0,1):
        for nu in(-1,0,1):
            if mu+nu!=kap: continue
            cg=CG(1,1,k,mu,nu,kap)
            if cg==0: continue
            for (mlk,msk,ck) in coupled(l,s,J,M):
                # L_1nu first (orbital), then J_1mu
                if lp!=l: continue
                mli=mlk+nu; msi=msk
                aL=L_op(l,mli,nu,l,mlk)
                if aL==0: continue
                for (mlf,msf,ampJ) in J_action(l,mli,msi,mu,s):
                    for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                        if mlb==mlf and msb==msf:
                            tot+=cb*ck*cg*aL*ampJ
    return simplify(tot)
def LJ_book(lp,s,Jp,Mp,k,kap,l,J,M):
    return simplify((-1)**(s+l-Jp+k+1)*CG(J,k,Jp,M,kap,Mp)*(2*J+1)
        *sqrt((2*k+1)*(2*l+1)*l*(l+1)*J*(J+1))*wigner_6j(1,1,k,Jp,J,J)*wigner_6j(l,s,J,Jp,1,l))
def JL_book(lp,s,Jp,Mp,k,kap,l,J,M):
    return simplify((-1)**(s+l-Jp+k+1)*CG(J,k,Jp,M,kap,Mp)
        *sqrt((2*k+1)*l*(l+1)*(2*l+1)*Jp*(Jp+1)*(2*Jp+1)*(2*J+1))
        *wigner_6j(1,1,k,Jp,J,Jp)*wigner_6j(l,s,J,Jp,1,l))
ccases=[(1,h,h,h,1,0,1,h,h),(1,h,3*h,h,1,0,1,3*h,h),(1,h,3*h,3*h,2,1,1,3*h,h),
        (2,h,5*h,h,1,0,2,3*h,h),(2,h,3*h,h,2,0,2,3*h,h),(1,h,3*h,h,2,-1,1,3*h,3*h)]
print("  86:",all(simplify(LJ_true(*c)-LJ_book(*c))==0 for c in ccases))
print("  87:",all(simplify(JL_true(*c)-JL_book(*c))==0 for c in ccases))
