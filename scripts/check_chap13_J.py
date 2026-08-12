"""Verify J matrix-element product formulas, Chap13 sec 13.2.5."""
from sympy.physics.wigner import clebsch_gordan as CG, gaunt, wigner_6j
from sympy import sqrt, Rational, simplify, S, pi

def Jsph(Mp,mu,J,M):           # <J Mp|J_{1mu}|J M>
    if mu==0:  return S(M) if Mp==M else S(0)
    if mu==1:  return -1/sqrt(2)*sqrt((J-M)*(J+M+1)) if Mp==M+1 else S(0)
    return  1/sqrt(2)*sqrt((J+M)*(J-M+1)) if Mp==M-1 else S(0)

print("=== eq 13.2.47: <J M'|{J1 x J1}_k,kappa|J M> ===")
def JJ_true(Mp,k,kap,J,M):
    t=S(0)
    for mu in(-1,0,1):
        nu=kap-mu
        if abs(nu)>1: continue
        t+=CG(1,1,k,mu,nu,kap)*sum(Jsph(Mp,mu,J,Mpp)*Jsph(Mpp,nu,J,M) for Mpp in range(-J,J+1))
    return simplify(t)
def JJ_book(Mp,k,kap,J,M):
    return simplify((-1)**(2*J+k)*J*(J+1)*sqrt((2*k+1)*(2*J+1))*wigner_6j(1,1,k,J,J,J)*CG(J,k,J,M,kap,Mp))
ok=all(simplify(JJ_true(Mp,k,kap,J,M)-JJ_book(Mp,k,kap,J,M))==0
       for J in(1,2) for k in(0,1,2) for M in range(-J,J+1)
       for kap in(-k,0,k) for Mp in[M+kap] if abs(Mp)<=J)
print("  match:",ok)

# --- coupled (l s J) basis, n-hat orbital rank-1, J total ---
def n_red(lp,Jp,l,J,s):        # eq 13.2.12  <l' s J'||n_1||l s J>  (s'=s)
    return (-1)**(s+J+lp+1)*sqrt((2*J+1)*(2*Jp+1)*(2*l+1))*CG(l,1,lp,0,0,0)*wigner_6j(l,s,J,Jp,1,lp)
def n_c(lp,Jp,Mp,mu,l,J,M,s):
    if lp<0 or Jp<0 or abs(Mp)>Jp: return S(0)
    return n_red(lp,Jp,l,J,s)/sqrt(2*Jp+1)*CG(J,1,Jp,M,mu,Mp)
def J_c(Jp,Mp,nu,J,M):         # J diagonal in l,s,J
    return sqrt(J*(J+1))*CG(J,1,J,M,nu,Mp) if Jp==J and abs(Mp)<=J else S(0)

def nJ_true(lp,Jp,Mp,k,kap,l,J,M,s):     # {n x J}_k
    t=S(0)
    for mu in(-1,0,1):
        nu=kap-mu
        if abs(nu)>1: continue
        t+=CG(1,1,k,mu,nu,kap)*n_c(lp,Jp,Mp,mu,l,J,M+nu,s)*J_c(J,M+nu,nu,J,M)
    return simplify(t)
def Jn_true(lp,Jp,Mp,k,kap,l,J,M,s):     # {J x n}_k
    t=S(0)
    for mu in(-1,0,1):
        nu=kap-mu
        if abs(nu)>1: continue
        t+=CG(1,1,k,mu,nu,kap)*J_c(Jp,Mp,mu,Jp,M+nu)*n_c(lp,Jp,M+nu,nu,l,J,M,s)
    return simplify(t)
def nJ_book(lp,Jp,Mp,k,kap,l,J,M,s):     # eq 13.2.52  (6j: J' J J)
    return simplify((-1)**(s+lp-Jp+k+1)*sqrt((2*l+1)*(2*k+1)*J*(J+1))*(2*J+1)
        *CG(l,1,lp,0,0,0)*wigner_6j(1,1,k,Jp,J,J)*wigner_6j(l,s,J,Jp,1,lp)*CG(J,k,Jp,M,kap,Mp))
def Jn_book(lp,Jp,Mp,k,kap,l,J,M,s):     # eq 13.2.56  (6j: J' J J' ; uses J'(J'+1))
    return simplify((-1)**(s+lp-Jp+k+1)*sqrt((2*l+1)*(2*k+1)*(2*Jp+1)*(2*J+1))*sqrt(Jp*(Jp+1))
        *CG(l,1,lp,0,0,0)*wigner_6j(1,1,k,Jp,J,Jp)*wigner_6j(l,s,J,Jp,1,lp)*CG(J,k,Jp,M,kap,Mp))

print("=== eq 13.2.52 {n x J} vs direct (6j = J' J J) ===")
h=Rational(1,2)
# (lp, Jp, Mp, k, kap, l, J, M, s=1/2), valid coupled states J=l+-1/2, l'=l+-1
cases=[(2,5*h,h,1,0,1,3*h,h,h),(0,h,h,1,0,1,3*h,h,h),(1,3*h,3*h,2,1,2,5*h,h,h),
       (2,3*h,h,1,0,1,h,h,h),(2,5*h,3*h,2,1,1,3*h,h,h)]
print("  ",all(simplify(nJ_true(*c)-nJ_book(*c))==0 for c in cases))
print("=== eq 13.2.56 {J x n} vs direct (6j = J' J J') ===")
print("  ",all(simplify(Jn_true(*c)-Jn_book(*c))==0 for c in cases))
