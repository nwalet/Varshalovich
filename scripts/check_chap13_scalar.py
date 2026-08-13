"""Verify scalar/vector products (Chap13 sec 13.2.9, eqs 13.2.150-174)."""
from sympy.physics.wigner import gaunt
from sympy import sqrt, Rational, S, pi, N, I
from sympy.physics.wigner import clebsch_gordan as _CG0, wigner_6j as _6j0
def CG(*a):
    try: return _CG0(*a)
    except Exception: return S(0)
def w6(*a):
    try: return _6j0(*a)
    except Exception: return S(0)
h=Rational(1,2)
def n_op(lp,mlp,mu,l,ml):
    if lp<0: return S(0)
    return sqrt(4*pi/3)*(-1)**mlp*gaunt(lp,1,l,-mlp,mu,ml)
def L_op(lp,mlp,mu,l,ml):
    if lp!=l or mlp!=ml+mu: return S(0)
    return sqrt(l*(l+1))*CG(l,1,l,ml,mu,mlp)
def S_op(sp,msp,mu,s,ms):
    if sp!=s or msp!=ms+mu: return S(0)
    return sqrt(s*(s+1))*CG(s,1,s,ms,mu,msp)
def gO(lp,mp,mu,l,m):
    t=S(0)
    for mu1 in (-1,0,1):
        mu2=mu-mu1
        if abs(mu2)<=1: t+=CG(1,1,1,mu1,mu2,mu)*n_op(lp,mp,mu1,l,m+mu2)*L_op(l,m+mu2,mu2,l,m)
    return -sqrt(2)*t
def coupled(l,s,J,M):
    out=[]
    for i in range(2*l+1):
        ml=-l+i; ms=M-ml
        if abs(ms)<=s:
            c=CG(l,s,J,ml,ms,M)
            if c!=0: out.append((ml,ms,c))
    return out
def close(a,b): return abs(complex(N(a-b)))<1e-9

# apply operator component to uncoupled |l ml, s ms> -> list (lf,mlf,msf,amp)
def apply(op,comp,l,ml,ms,s):
    r=[]
    if op=='n':
        for lf in (l-1,l+1):
            a=n_op(lf,ml+comp,comp,l,ml)
            if a!=0: r.append((lf,ml+comp,ms,a))
    elif op=='g':
        for lf in (l-1,l+1):
            a=gO(lf,ml+comp,comp,l,ml)
            if a!=0: r.append((lf,ml+comp,ms,a))
    elif op=='L':
        a=L_op(l,ml+comp,comp,l,ml)
        if a!=0: r.append((l,ml+comp,ms,a))
    elif op=='S':
        a=S_op(s,ms+comp,comp,s,ms)
        if a!=0: r.append((l,ml,ms+comp,a))
    elif op=='J':
        a=L_op(l,ml+comp,comp,l,ml)
        if a!=0: r.append((l,ml+comp,ms,a))
        b=S_op(s,ms+comp,comp,s,ms)
        if b!=0: r.append((l,ml,ms+comp,b))
    return r

# <l' s J' M'| A_1al B_1be |l s J M>  (B acts first)
def AB(A,al,B,be,lp,s,Jp,Mp,l,J,M):
    tot=S(0)
    for (mlk,msk,ck) in coupled(l,s,J,M):
        for (li,mli,msi,ampB) in apply(B,be,l,mlk,msk,s):
            for (lf,mlf,msf,ampA) in apply(A,al,li,mli,msi,s):
                if lf!=lp: continue
                for (mlb,msb,cb) in coupled(lp,s,Jp,Mp):
                    if mlb==mlf and msb==msf: tot+=cb*ck*ampB*ampA
    return tot
def scalar(A,B,lp,s,Jp,Mp,l,J,M):
    return sum((-1)**mu*AB(A,mu,B,-mu,lp,s,Jp,Mp,l,J,M) for mu in(-1,0,1))
def vec(A,B,kap,lp,s,Jp,Mp,l,J,M):   # [AxB]_kap = -i sqrt2 {A x B}_{1kap}
    t=S(0)
    for al in(-1,0,1):
        be=kap-al
        if abs(be)>1: continue
        t+=CG(1,1,1,al,be,kap)*AB(A,al,B,be,lp,s,Jp,Mp,l,J,M)
    return -I*sqrt(2)*t

# orbital-only wrappers (s=0 dummy, J=l, M=m)
def scal_orb(A,B,lp,l,m): return scalar(A,B,lp,S(0),lp,m,l,l,m)
def vec_orb(A,B,kap,lp,l,m): return vec(A,B,kap,lp,S(0),lp,m+kap,l,l,m)

print("=== check convention: [n x grad]_kappa = i L_kappa  (eq 13.2.151) ===")
ok=all(close(vec_orb('n','g',kap,lp,l,m), I*L_op(lp,m+kap,kap,l,m))
       for (l,lp,m,kap) in [(1,1,0,1),(2,2,1,0),(1,1,-1,1),(2,2,0,-1)])
print("  ",ok)

# ---- 13.2.152 (n.J) ----
def nJ_book(lp,s,Jp,Mp,l,J,M):
    if Jp!=J or Mp!=M: return S(0)
    return (-1)**(s+lp+J+1)*sqrt(J*(J+1)*(2*J+1)*(2*l+1))*CG(l,1,lp,0,0,0)*w6(l,s,J,J,1,lp)
c=[(1,h,3*h,h,2,3*h,h),(2,h,3*h,h,1,3*h,h),(1,h,h,h,2,h,h),(3,h,5*h,h,2,5*h,h)]
print("152 (n.J):",all(close(scalar('n','J',*x),nJ_book(*x)) for x in c))

# ---- 13.2.153 [n x J] ----
def nxJ_book(lp,s,Jp,Mp,kap,l,J,M):
    return (I*h*(-1)**(J+s+lp)*(Jp*(Jp+1)-J*(J+1)-2)*sqrt((2*l+1)*(2*J+1))
        *CG(l,1,lp,0,0,0)*w6(l,s,J,Jp,1,lp)*CG(J,1,Jp,M,kap,Mp))
cv=[(2,h,5*h,h,0,1,3*h,h),(0,h,h,h,0,1,3*h,h),(2,h,3*h,h,0,1,h,h),(2,h,3*h,3*h,1,1,h,-h)]
print("153 [n x J]:",all(close(vec('n','J',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),nxJ_book(*x)) for x in cv))

# ---- 13.2.154 [J x n]: test C^{l'0}_{l0,10} (first arg l) vs 10 ----
def Jxn_book(lp,s,Jp,Mp,kap,l,J,M,first):
    return (I*h*(-1)**(J+s+lp+1)*(Jp*(Jp+1)-J*(J+1)+2)*sqrt((2*l+1)*(2*J+1))
        *CG(first,1,lp,0,0,0)*w6(l,s,J,Jp,1,lp)*CG(J,1,Jp,M,kap,Mp))
print("154 [J x n] first=l:",all(close(vec('J','n',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),Jxn_book(*x,x[5])) for x in cv))
print("154 [J x n] first=1:",all(close(vec('J','n',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),Jxn_book(*x,1)) for x in cv))

# ---- 13.2.164/165 (L x J), (J x L): test 6j last arg l vs l' ----
def LxJ_book(lp,s,Jp,Mp,kap,l,J,M,last):
    return (I*h*(-1)**(s+l+J+1)*(J*(J+1)-Jp*(Jp+1)+2)*sqrt((2*J+1)*l*(l+1)*(2*l+1))
        *w6(l,s,J,Jp,1,last)*CG(J,1,Jp,M,kap,Mp)) if lp==l else S(0)
cL=[(1,h,3*h,h,0,1,3*h,h),(1,h,h,h,0,1,3*h,h),(2,h,5*h,h,0,2,3*h,h),(1,h,3*h,3*h,1,1,h,-h)]
print("164 [L x J] 6j last=l :",all(close(vec('L','J',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),LxJ_book(*x,x[5])) for x in cL))
def JxL_book(lp,s,Jp,Mp,kap,l,J,M,last):
    return (I*h*(-1)**(s+l+J+1)*(Jp*(Jp+1)-J*(J+1)+2)*sqrt((2*J+1)*l*(l+1)*(2*l+1))
        *w6(l,s,J,Jp,1,last)*CG(J,1,Jp,M,kap,Mp)) if lp==l else S(0)
print("165 [J x L] 6j last=l :",all(close(vec('J','L',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),JxL_book(*x,x[5])) for x in cL))
print("165 [J x L] 6j last=l':",all(close(vec('J','L',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),JxL_book(*x,x[0])) for x in cL))

# ---- 13.2.166 (n.S) 6j {1 s s; J l l'} ----
def nS_book(lp,s,Jp,Mp,l,J,M):
    if Jp!=J or Mp!=M: return S(0)
    return (-1)**(l+s+J)*sqrt(s*(s+1)*(2*s+1)*(2*l+1))*CG(l,1,lp,0,0,0)*w6(1,s,s,J,l,lp)
print("166 (n.S):",all(close(scalar('n','S',*x),nS_book(*x)) for x in c))

# ---- 13.2.169 [n x S]: 6j {J' l' 3->s; l J 1} ----
def nxS_book(lp,s,Jp,Mp,kap,l,J,M,mid):
    return (I*(-1)**(lp+J+s)*sqrt((2*l+1)*(2*J+1))*CG(l,1,lp,0,0,0)
        *h*((Jp-lp)*(Jp+lp+1)-(J-l)*(J+l+1))*w6(Jp,lp,mid,l,J,1)*CG(J,1,Jp,M,kap,Mp))
print("169 [n x S] 6j mid=s:",all(close(vec('n','S',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),nxS_book(*x,x[1])) for x in cv))
print("169 [n x S] 6j mid=3:",all(close(vec('n','S',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),nxS_book(*x,3)) for x in cv))

# ---- 13.2.173/174 [S x J],[J x S] 6j {s l J; J' 1 s} ----
def SxJ_book(lp,s,Jp,Mp,kap,l,J,M):
    return (I*h*(-1)**(l+s+Jp+1)*(J*(J+1)-Jp*(Jp+1)+2)*sqrt((2*J+1)*s*(s+1)*(2*s+1))
        *w6(s,l,J,Jp,1,s)*CG(J,1,Jp,M,kap,Mp)) if lp==l else S(0)
print("173 [S x J]:",all(close(vec('S','J',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),SxJ_book(*x)) for x in cL))
def JxS_book(lp,s,Jp,Mp,kap,l,J,M):
    return (I*h*(-1)**(l+s+Jp+1)*(Jp*(Jp+1)-J*(J+1)+2)*sqrt((2*J+1)*s*(s+1)*(2*s+1))
        *w6(s,l,J,Jp,1,s)*CG(J,1,Jp,M,kap,Mp)) if lp==l else S(0)
print("174 [J x S]:",all(close(vec('J','S',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),JxS_book(*x)) for x in cL))

# ============ nabla_Omega products: resolve bracket [coef*sqrt((l+1)/(2l+3)) d_{l',l+1} + (l+1)sqrt(l/(2l-1)) d_{l',l-1}] ============
def brk(coef1,lp,l):
    if lp==l+1: return coef1*sqrt(S(l+1)/(2*l+3))
    if lp==l-1: return (l+1)*sqrt(S(l)/(2*l-1))
    return S(0)
# 154 discriminating (l=2)
cv2=[(3,h,5*h,h,0,2,3*h,h),(1,h,h,h,0,2,3*h,h),(3,h,7*h,h,0,2,5*h,h)]
def Jxn_b(lp,s,Jp,Mp,kap,l,J,M,first):
    return (I*h*(-1)**(J+s+lp+1)*(Jp*(Jp+1)-J*(J+1)+2)*sqrt((2*l+1)*(2*J+1))
        *CG(first,1,lp,0,0,0)*w6(l,s,J,Jp,1,lp)*CG(J,1,Jp,M,kap,Mp))
print("154 disc first=l:",all(close(vec('J','n',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),Jxn_b(*x,x[5])) for x in cv2))
print("154 disc first=1:",all(close(vec('J','n',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),Jxn_b(*x,1)) for x in cv2))

# 13.2.155 (grad.J)
def gJ_book(lp,s,Jp,Mp,l,J,M,coef1):
    if Jp!=J or Mp!=M: return S(0)
    return (-1)**(s+lp+J)*sqrt((2*lp+1)*(2*J+1)*J*(J+1))*brk(coef1,lp,l)*w6(l,s,J,J,1,lp)
cg155=[(2,h,3*h,h,1,3*h,h),(0,h,h,h,1,h,h),(3,h,5*h,h,2,5*h,h),(1,h,h,h,2,3*h,h)]
print("155 (grad.J) coef1=l:",all(close(scalar('g','J',*x),gJ_book(*x,x[4]-1 if False else x[4])) for x in cg155))  # placeholder
print("155 (grad.J) coef1=l:",all(close(scalar('g','J',*x),gJ_book(*x,x[4])) for x in cg155))
print("155 (grad.J) coef1=1:",all(close(scalar('g','J',*x),gJ_book(*x,1)) for x in cg155))

# 13.2.156 [grad x J]
def gxJ_book(lp,s,Jp,Mp,kap,l,J,M,coef1):
    return (I*h*(-1)**(J+s+lp+1)*sqrt((2*lp+1)*(2*J+1))*(Jp*(Jp+1)-J*(J+1)-2)
        *brk(coef1,lp,l)*w6(l,s,J,Jp,1,lp)*CG(J,1,Jp,M,kap,Mp))
cvg=[(2,h,5*h,h,0,1,3*h,h),(0,h,h,h,0,1,3*h,h),(2,h,3*h,h,0,1,h,h),(3,h,5*h,h,0,2,3*h,h)]
print("156 [grad x J] coef1=l:",all(close(vec('g','J',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),gxJ_book(*x,x[5])) for x in cvg))
print("156 [grad x J] coef1=1:",all(close(vec('g','J',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),gxJ_book(*x,1)) for x in cvg))

# 13.2.161 [grad x L] (orbital)
def gxL_book(lp,l,m,kap,coef1):
    return (I*h*(lp*(lp+1)-l*(l+1)-2)*brk(coef1,lp,l)*CG(l,1,lp,m,kap,m+kap))
cgL=[(2,1,0,1),(0,1,0,1),(3,2,1,-1),(1,2,0,1)]
print("161 [grad x L] coef1=1:",all(close(vec_orb('g','L',kap,lp,l,m),gxL_book(lp,l,m,kap,1)) for (lp,l,m,kap) in cgL))
print("161 [grad x L] coef1=l:",all(close(vec_orb('g','L',kap,lp,l,m),gxL_book(lp,l,m,kap,l)) for (lp,l,m,kap) in cgL))

# 13.2.167 (grad.S)
def gS_book(lp,s,Jp,Mp,l,J,M,coef1):
    if Jp!=J or Mp!=M: return S(0)
    return (-1)**(l+s+J+1)*sqrt(s*(s+1)*(2*s+1)*(2*lp+1))*brk(coef1,lp,l)*w6(1,s,s,J,l,lp)
print("167 (grad.S) coef1=1:",all(close(scalar('g','S',*x),gS_book(*x,1)) for x in cg155))
print("167 (grad.S) coef1=l:",all(close(scalar('g','S',*x),gS_book(*x,x[4])) for x in cg155))

# 13.2.157 [J x grad] : phase ambiguous (garbled J++l')
def Jxg_book(lp,s,Jp,Mp,kap,l,J,M,ph):
    return (I*h*(-1)**ph*(Jp*(Jp+1)-J*(J+1)+2)*sqrt((2*lp+1)*(2*J+1))
        *brk(l,lp,l)*w6(l,s,J,Jp,1,lp)*CG(J,1,Jp,M,kap,Mp))
for nm,f in [("J+s+l'",lambda x:x[6]+x[1]+x[0]),("J+s+l'+1",lambda x:x[6]+x[1]+x[0]+1)]:
    print(f"157 [J x grad] phase {nm}:",all(close(vec('J','g',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),Jxg_book(*x,f(x))) for x in cvg))

# 13.2.162 [L x grad] (orbital) : coef1 and phase sign
def Lxg_book(lp,l,m,kap,coef1,sgn):
    return (sgn*I*h*(lp*(lp+1)-l*(l+1)+2)*brk(coef1,lp,l)*CG(l,1,lp,m,kap,m+kap))
print("162 [L x grad] coef1=l, sign=-:",all(close(vec_orb('L','g',kap,lp,l,m),Lxg_book(lp,l,m,kap,l,-1)) for (lp,l,m,kap) in cgL))

# 13.2.170 [grad x S] : phase iota->s
def gxS_book(lp,s,Jp,Mp,kap,l,J,M,ph):
    return (I*(-1)**ph*sqrt((2*lp+1)*(2*J+1))*w6(Jp,lp,s,l,J,1)
        *h*((Jp-lp)*(Jp+lp+1)-(J-l)*(J+l+1))*brk(l,lp,l)*CG(J,1,Jp,M,kap,Mp))
print("170 [grad x S] phase l'+J+s+1:",all(close(vec('g','S',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),gxS_book(*x,x[0]+x[6]+x[1]+1)) for x in cvg))

# 13.2.171 [L x S] (orbital-diagonal, coupled)
def LxS_book(lp,s,Jp,Mp,kap,l,J,M):
    return (I*(-1)**(l+J+s)*sqrt((2*J+1)*l*(l+1)*(2*l+1))
        *h*((Jp-l)*(Jp+l+1)-(J-l)*(J+l+1))*w6(Jp,l,s,l,J,1)*CG(J,1,Jp,M,kap,Mp)) if lp==l else S(0)
print("171 [L x S]:",all(close(vec('L','S',x[4],x[0],x[1],x[2],x[3],x[5],x[6],x[7]),LxS_book(*x)) for x in cL))
