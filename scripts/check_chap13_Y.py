"""Verify spherical-harmonic operator Y_{L nu} matrix elements, Chap13 sec 13.2.8 (part A: 104-127)."""
from sympy.physics.wigner import clebsch_gordan as CG, gaunt, wigner_6j, wigner_9j
from sympy import sqrt, Rational, simplify, S, pi

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
def gO(lp,mp,mu,l,m):                 # grad_Omega = -i[n x L]
    t=S(0)
    for mu1 in (-1,0,1):
        mu2=mu-mu1
        if abs(mu2)<=1: t+=CG(1,1,1,mu1,mu2,mu)*n_op(lp,mp,mu1,l,m+mu2)*L_op(l,m+mu2,mu2,l,m)
    return simplify(-sqrt(2)*t)

# ---------- 13.2.108  n_1mu Y_Lnu  (direct, orbital) ----------
print("=== 13.2.108  n_1mu Y_Lnu ===")
def nY_true(lp,mp,mu,L,nu,l,m):
    return simplify(sum(n_op(lp,mp,mu,lpp,m+nu)*Y_op(lpp,m+nu,L,nu,l,m) for lpp in (lp-1,lp+1)))
def nY_book(lp,mp,mu,L,nu,l,m):
    pre=-sqrt((2*L+1)*(2*l+1)/(4*pi*(2*lp+1)))
    t1=sqrt(S(lp+1)/(2*lp+3))*CG(l,L,lp+1,0,0,0)*CG(l,L,lp+1,m,nu,m+nu)*CG(lp+1,1,lp,m+nu,mu,mp)
    t2=sqrt(S(lp)/(2*lp-1))*CG(l,L,lp-1,0,0,0)*CG(l,L,lp-1,m,nu,m+nu)*CG(lp-1,1,lp,m+nu,mu,mp)
    return simplify(pre*(t1-t2))
cn=[(2,1,1,2,0,1,0),(1,0,1,2,1,2,0),(3,2,-1,1,0,2,1),(2,0,1,1,-1,2,1),(1,1,0,2,0,1,1)]
print("  ",all(simplify(nY_true(*c)-nY_book(*c))==0 for c in cn))

# ---------- 13.2.109  {n x Y_L}_{L' nu'} ----------
print("=== 13.2.109  {n x Y_L}_{L'nu'} ===")
def nYt_true(lp,mp,Lp,vp,L,l,m):      # tensor product rank L'
    t=S(0)
    for mu in(-1,0,1):
        nu=vp-mu
        for lpp in (lp-1,lp+1):
            t+=CG(1,L,Lp,mu,nu,vp)*n_op(lp,mp,mu,lpp,m+nu)*Y_op(lpp,m+nu,L,nu,l,m)
    return simplify(t)
def nYt_book(lp,mp,Lp,vp,L,l,m):
    pre=(-1)**(Lp+lp+l)*sqrt((2*L+1)*(2*Lp+1)*(2*l+1)/(4*pi*(2*lp+1)))*CG(l,Lp,lp,m,vp,mp)
    s=S(0)
    for k in range(abs(l-L),l+L+1):
        s+=sqrt(2*k+1)*CG(k,1,lp,0,0,0)*CG(l,L,k,0,0,0)*wigner_6j(L,1,Lp,lp,l,k)
    return pre*s
from sympy import N as _N
def _close(a,b): return abs(complex(_N(a-b)))<1e-9
cnt=[(2,1,2,1,2,2,0),(1,0,1,0,1,1,0),(3,1,2,0,2,2,1),(2,-1,3,0,2,3,1)]
print("  ",all(_close(nYt_true(*c),nYt_book(*c)) for c in cnt))

# ---------- 13.2.110  grad_Omega Y_Lnu ----------
print("=== 13.2.110  (grad_O)_1mu Y_Lnu ===")
def gY_true(lp,mp,mu,L,nu,l,m):
    return simplify(sum(gO(lp,mp,mu,lpp,m+nu)*Y_op(lpp,m+nu,L,nu,l,m) for lpp in (lp-1,lp+1)))
def gY_book(lp,mp,mu,L,nu,l,m):
    pre=-sqrt((2*L+1)*(2*l+1)/(4*pi*(2*lp+1)))
    t1=(lp-1)*sqrt(S(lp)/(2*lp-1))*CG(l,L,lp-1,0,0,0)*CG(l,L,lp-1,m,nu,m+nu)*CG(lp-1,1,lp,m+nu,mu,mp)
    t2=(lp+2)*sqrt(S(lp+1)/(2*lp+3))*CG(l,L,lp+1,0,0,0)*CG(l,L,lp+1,m,nu,m+nu)*CG(lp+1,1,lp,m+nu,mu,mp)
    return simplify(pre*(t1+t2))
print("  ",all(simplify(gY_true(*c)-gY_book(*c))==0 for c in cn))

# ---------- 13.2.116  L_1mu Y_Lnu ----------
print("=== 13.2.116  L_1mu Y_Lnu ===")
def LY_true(lp,mp,mu,L,nu,l,m):
    return simplify(L_op(lp,mp,mu,lp,m+nu)*Y_op(lp,m+nu,L,nu,l,m))
def LY_book(lp,mp,mu,L,nu,l,m):
    return simplify(sqrt((2*L+1)*(2*l+1)*lp*(lp+1)/(4*pi*(2*lp+1)))
        *CG(l,L,lp,0,0,0)*CG(l,L,lp,m,nu,m+nu)*CG(lp,1,lp,m+nu,mu,mp))
print("  ",all(simplify(LY_true(*c)-LY_book(*c))==0 for c in cn))

# ---------- 13.2.117  {L x Y_L}_{L'nu'} ----------
print("=== 13.2.117  {L x Y_L}_{L'nu'} ===")
def LYt_true(lp,mp,Lp,vp,L,l,m):
    t=S(0)
    for mu in(-1,0,1):
        nu=vp-mu
        t+=CG(1,L,Lp,mu,nu,vp)*L_op(lp,mp,mu,lp,m+nu)*Y_op(lp,m+nu,L,nu,l,m)
    return simplify(t)
def LYt_book(lp,mp,Lp,vp,L,l,m):
    return simplify((-1)**(Lp+lp+l)*sqrt((2*Lp+1)*(2*L+1)*(2*l+1)*(lp+1)*lp/(4*pi))
        *CG(l,Lp,lp,m,vp,mp)*CG(l,L,lp,0,0,0)*wigner_6j(L,1,Lp,lp,l,lp))
print("  ",all(simplify(LYt_true(*c)-LYt_book(*c))==0 for c in cnt))

# ---------- 13.2.120 / 121  commutator [L_1mu, Y_Lnu] ----------
print("=== 13.2.118/120  R(L,Y) = sqrt(L(L+1)) C^{L,mu+nu}_{L nu,1 mu} Y_{L,mu+nu} ===")
def comm_true(lp,mp,mu,L,nu,l,m):     # <l'm'|[L_1mu,Y_Lnu]|lm>
    a=sum(L_op(lp,mp,mu,lp,mpp)*Y_op(lp,mpp,L,nu,l,m) for mpp in [m+nu])  # L then Y? careful
    # [L_1mu, Y_Lnu] = L_1mu Y_Lnu - Y_Lnu L_1mu
    lhs=L_op(lp,mp,mu,lp,m+nu)*Y_op(lp,m+nu,L,nu,l,m)
    rhs=sum(Y_op(lp,mp,L,nu,l+0,m+mu)*0 for _ in [0])  # placeholder
    # Y_Lnu L_1mu : L acts first on |lm> (diagonal), then Y
    rhs=Y_op(lp,mp,L,nu,l,m+mu)*L_op(l,m+mu,mu,l,m)
    return simplify(lhs-rhs)
def comm_book118(lp,mp,mu,L,nu,l,m):
    return simplify(sqrt(L*(L+1))*CG(L,1,L,nu,mu,mu+nu)*Y_op(lp,mp,L,mu+nu,l,m))
cc=[(2,1,1,2,0,1,0),(1,0,-1,1,1,2,1),(3,2,1,2,-1,3,2),(2,0,1,1,0,2,1)]
print("  118:",all(simplify(comm_true(*c)-comm_book118(*c))==0 for c in cc))
def comm120_book(lp,mp,mu,L,nu,l,m):
    return simplify(sqrt(L*(L+1)*(2*L+1)*(2*l+1)/(4*pi*(2*lp+1)))
        *CG(l,L,lp,m,mu+nu,mp)*CG(L,1,L,nu,mu,nu+mu)*CG(l,L,lp,0,0,0))
print("  120:",all(simplify(comm_true(*c)-comm120_book(*c))==0 for c in cc))

# ================= coupled (l s J) with Y (part A: 122,124) =================
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

print("=== 13.2.122  S_1mu Y_Lnu  (coupled) : test which CG ===")
def SY_true(lp,s,Jp,Mp,mu,L,nu,l,J,M):
    tot=S(0)
    for (mlk,msk,ck) in coupled(l,s,J,M):
        for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
            amp=Y_op(lp,mlb,L,nu,l,mlk)*S_op(s,msb,mu,s,msk)
            if amp!=0: tot+=cb*ck*amp
    return simplify(tot)
def SY_book(lp,s,Jp,Mp,mu,L,nu,l,J,M,cgfirst):
    pre=sqrt((2*J+1)*(2*l+1)*(2*L+1)*s*(s+1)*(2*s+1)/(4*pi))*CG(cgfirst,0,(l if cgfirst==1 else 1)*0+ (L if cgfirst==l else L),0,lp,0)
    # build CG explicitly below instead
    return None
def SY_book_l(lp,s,Jp,Mp,mu,L,nu,l,J,M):
    pre=sqrt((2*J+1)*(2*l+1)*(2*L+1)*s*(s+1)*(2*s+1)/(4*pi))*CG(l,L,lp,0,0,0)
    ss=S(0)
    for k in range(0, L+2):
        ss+=(-1)**(L+k+1)*sqrt(2*k+1)*CG(J,k,Jp,M,mu+nu,Mp)*CG(1,L,k,mu,nu,mu+nu)*wigner_9j(l,L,lp,s,1,s,J,k,Jp)
    return simplify(pre*ss)
def SY_book_1(lp,s,Jp,Mp,mu,L,nu,l,J,M):
    pre=sqrt((2*J+1)*(2*l+1)*(2*L+1)*s*(s+1)*(2*s+1)/(4*pi))*CG(1,0,L,0,lp,0)
    ss=S(0)
    for k in range(0, L+2):
        ss+=(-1)**(L+k+1)*sqrt(2*k+1)*CG(J,k,Jp,M,mu+nu,Mp)*CG(1,L,k,mu,nu,mu+nu)*wigner_9j(l,L,lp,s,1,s,J,k,Jp)
    return simplify(pre*ss)
csy=[(2,h,5*h,h,0,2,0,2,3*h,h),(1,h,3*h,h,1,1,0,1,h,h),(0,h,h,h,0,1,0,1,3*h,h),(2,h,3*h,3*h,1,1,1,1,h,h)]
print("  122 (CG first=l):",all(simplify(SY_true(*c)-SY_book_l(*c))==0 for c in csy))
print("  122 (CG first=1):",all(simplify(SY_true(*c)-SY_book_1(*c))==0 for c in csy))

print("=== 13.2.124  J_1mu Y_Lnu  (coupled) : test which CG ===")
def JY_true(lp,s,Jp,Mp,mu,L,nu,l,J,M):
    tot=S(0)
    for (mlk,msk,ck) in coupled(l,s,J,M):
        # Y first (orbital) -> intermediate (mli=mlk+nu, msk); then J_1mu
        mli=mlk+nu
        aY=Y_op(lp,mli,L,nu,l,mlk)
        if aY==0: continue
        for (mlf,msf,ampJ) in J_action(lp,mli,msk,mu,s):
            for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                if mlb==mlf and msb==msf: tot+=cb*ck*aY*ampJ
    return simplify(tot)
def JY_book(lp,s,Jp,Mp,mu,L,nu,l,J,M,cg):
    return simplify((-1)**(lp+s+J+L)*sqrt((2*L+1)*(2*l+1)*(2*J+1)*Jp*(Jp+1)/(4*pi))
        *CG(J,L,Jp,M,nu,M+nu)*CG(Jp,1,Jp,M+nu,mu,Mp)*cg*wigner_6j(l,s,J,Jp,L,lp))
cjy=[(2,h,5*h,h,0,2,0,2,3*h,h),(1,h,3*h,h,1,1,0,1,h,h),(2,h,3*h,3*h,1,1,1,2,h,h),(0,h,h,h,0,1,0,1,h,h)]
okl=all(simplify(JY_true(*c)-JY_book(*c,CG(l,L,lp,0,0,0)))==0 for (lp,s,Jp,Mp,mu,L,nu,l,J,M) in cjy for c in [(lp,s,Jp,Mp,mu,L,nu,l,J,M)])
ok1=all(simplify(JY_true(*c)-JY_book(*c,CG(L,1,lp,0,0,0)))==0 for (lp,s,Jp,Mp,mu,L,nu,l,J,M) in cjy for c in [(lp,s,Jp,Mp,mu,L,nu,l,J,M)])
print("  124 (CG=C^{l'0}_{l0,L0}):",okl)
print("  124 (CG=C^{l'0}_{L0,10}):",ok1)
