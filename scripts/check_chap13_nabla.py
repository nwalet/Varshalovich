"""Verify grad_Omega matrix elements, Chap13 sec 13.2.4.
Ground truth: grad_Omega = -i[n x L] built from n-hat (Gaunt) and L.
(grad_O)_1mu = -sqrt2 sum_{mu1 mu2} C^{1mu}_{1mu1,1mu2} n_1mu1 L_1mu2."""
from sympy.physics.wigner import clebsch_gordan as CG, gaunt
from sympy import sqrt, Rational, simplify, S, pi

def n_sph(lp,mp,mu,l,m):
    return S(0) if lp<0 else sqrt(4*pi/3)*(-1)**mp*gaunt(lp,1,l,-mp,mu,m)
def L_sph(mp,mu,l,m):
    if mu==0:  return S(m) if mp==m else S(0)
    if mu==1:  return -1/sqrt(2)*sqrt((l-m)*(l+m+1)) if mp==m+1 else S(0)
    return  1/sqrt(2)*sqrt((l+m)*(l-m+1)) if mp==m-1 else S(0)
def gO(lp,mp,mu,l,m):
    t=S(0)
    for mu1 in (-1,0,1):
        mu2=mu-mu1
        if abs(mu2)<=1: t+=CG(1,1,1,mu1,mu2,mu)*n_sph(lp,mp,mu1,l,m+mu2)*L_sph(m+mu2,mu2,l,m)
    return simplify(-sqrt(2)*t)

def red(lp,l):                       # eq 13.2.27
    if lp==l+1: return -l*sqrt(l+1)
    if lp==l-1: return -(l+1)*sqrt(S(l))
    return S(0)
def gO_book(lp,mp,mu,l,m):           # eq 13.2.26/27
    return S(0) if lp<0 else red(lp,l)/sqrt(2*lp+1)*CG(l,1,lp,m,mu,mp)

print("=== eq 13.2.26/27 reduced ME vs independent grad_O = -i[n x L] ===")
ok=all(simplify(gO(lp,mp,mu,l,m)-gO_book(lp,mp,mu,l,m))==0
       for (lp,mp,mu,l,m) in [(2,1,1,1,0),(0,0,0,1,0),(3,1,-1,2,2),(1,0,0,2,0),(2,0,-1,3,1)])
print("  match:", ok)

print("=== (grad_O . grad_O) = -l(l+1) ===")
def scalar(l,m):
    return simplify(sum((-1)**mu*gO(l,m,mu,lpp,m-mu)*gO(lpp,m-mu,-mu,l,m)
                        for mu in (-1,0,1) for lpp in (l-1,l,l+1)))
print("  ", all(scalar(l,m)==-l*(l+1) for l in (1,2,3) for m in range(-l,l+1)))

print("=== eq 13.2.29 scalar + rank-2 terms (book) vs true, on the pure-rank-2")
print("    channels kappa=+-2 and l'=l+-2 (where no rank-1 contributes) ===")
def prod_true(lp,mp,mu,nu,l,m):
    return simplify(sum(gO(lp,mp,mu,lpp,m+nu)*gO(lpp,m+nu,nu,l,m) for lpp in (l-1,l,l+1)))
def book_r02(lp,mp,mu,nu,l,m):       # eq 13.2.29 as printed (rank-0 scalar + rank-2)
    k=mu+nu
    if mp!=m+k: return S(0)
    out=S(0)
    if lp==l and mu==-nu and mp==m: out+=Rational((-1)**(1-mu),3)*l*(l+1)
    pref=CG(1,1,2,mu,nu,k)
    if lp==l+2: out+=pref*(l*(l+1)*sqrt((l+1)*(l+2))/sqrt((2*l+3)*(2*l+5)))*CG(l,2,l+2,m,k,mp)
    if lp==l:   out+=pref*(Rational(1,2*l+1)*sqrt(Rational(l*(l+1),6*(2*l-1)*(2*l+3)))*(4*l**3+6*l**2-4*l-3))*CG(l,2,l,m,k,mp)
    if lp==l-2: out+=pref*(l*(l+1)*sqrt(l*(l-1))/sqrt((2*l-1)*(2*l-3)))*CG(l,2,l-2,m,k,mp)
    return simplify(out)
pure=[(4,2,1,1,2,0),(2,2,1,1,2,0),(3,3,1,1,3,1),(1,0,1,-1,3,0),(2,-2,-1,-1,2,0),(2,0,1,-1,2,1)]
print("  ", all(simplify(prod_true(*c)-book_r02(*c))==0 for c in pure))

print("=== NOTE: on the kappa=+-1, l'=l diagonal channel the true product carries")
print("    an extra rank-1 (commutator) piece absent from 13.2.29: ===")
for (l,m) in [(2,0),(3,0)]:
    extra=simplify(prod_true(l,m+1,1,0,l,m)-book_r02(l,m+1,1,0,l,m))
    print(f"    l={l} m={m}: true - book(13.2.29) = {extra}   (nonzero => rank-1 term omitted)")
