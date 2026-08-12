"""Verify n-hat (unit vector) matrix elements, Chap13 sec 13.2.3."""
from sympy.physics.wigner import clebsch_gordan, gaunt
from sympy import sqrt, simplify, pi, Rational, nsimplify

def me_exact(lp, mp, mu, l, m):
    """<l' m'| n_{1 mu} | l m>, with n_{1mu} = sqrt(4pi/3) Y_{1mu}.
    <l'm'|Y_{1mu}|lm> = (-1)^{m'} * gaunt(l',1,l,-m',mu,m)."""
    return sqrt(4*pi/3) * (-1)**mp * gaunt(lp, 1, l, -mp, mu, m)

def me_formula(lp, mp, mu, l, m, first):
    """eq 13.2.10 with reduced ME 13.2.11 = sqrt(2l+1) C^{l'0}_{first 0,1 0}."""
    red = sqrt(2*l+1) * clebsch_gordan(first, 1, lp, 0, 0, 0)
    return red/sqrt(2*lp+1) * clebsch_gordan(l, 1, lp, m, mu, mp)

print("=== eq 13.2.10/11: which CG first-arg reproduces <l'm'|n_1mu|lm>? ===")
cases = [(2,1,1,1,0),(3,2,1,2,1),(2,0,0,3,0),(4,-1,-1,3,0),(1,1,1,2,0),(3,1,-1,2,2)]
for (lp,mp,mu,l,m) in cases:
    ex = simplify(me_exact(lp,mp,mu,l,m))
    f_l = simplify(me_formula(lp,mp,mu,l,m,l))   # C^{l'0}_{l0,10}  (proposed fix)
    f_1 = simplify(me_formula(lp,mp,mu,l,m,1))   # C^{l'0}_{10,10}  (current tex)
    print(f"  l'={lp} m'={mp} mu={mu:+d} l={l} m={m}: exact={ex}  first=l? {ex==f_l}  first=1? {ex==f_1}")

print()
print("=== eq 13.2.13: component matrix elements (take upper sign) ===")
from sympy import symbols
# line 1: <l+1, m+1| n_{1,+1}|l m> = sqrt((l+m+1)(l+m+2)/(2(2l+1)(2l+3)))
def num(l,m):
    return simplify(me_exact(l+1, m+1, 1, l, m))
def book13a(l,m):
    return sqrt(Rational((l+m+1)*(l+m+2),2*(2*l+1)*(2*l+3)))
def book13b(l,m):  # <l+1,m|n_{10}|lm> = sqrt((l-m+1)(l+m+1)/((2l+1)(2l+3)))
    return sqrt(Rational((l-m+1)*(l+m+1),(2*l+1)*(2*l+3)))
def num10(l,m):
    return simplify(me_exact(l+1, m, 0, l, m))
def book13c(l,m):  # <l-1,m+1|n_{1,+1}|lm> = -sqrt((l-m-1)(l-m)/(2(2l+1)(2l-1)))
    return -sqrt(Rational((l-m-1)*(l-m),2*(2*l+1)*(2*l-1)))
def numm1(l,m):
    return simplify(me_exact(l-1, m+1, 1, l, m))
def book13d(l,m):  # <l-1,m|n_{10}|lm> = sqrt((l-m)(l+m)/((2l+1)(2l-1)))
    return sqrt(Rational((l-m)*(l+m),(2*l+1)*(2*l-1)))
def numm10(l,m):
    return simplify(me_exact(l-1, m, 0, l, m))
for (l,m) in [(2,1),(3,0),(3,-2),(4,2)]:
    print(f"  l={l} m={m}: "
          f"13a {simplify(num(l,m)-book13a(l,m))==0}  "
          f"13b {simplify(num10(l,m)-book13b(l,m))==0}  "
          f"13c {simplify(numm1(l,m)-book13c(l,m))==0}  "
          f"13d {simplify(numm10(l,m)-book13d(l,m))==0}")

print()
print("=== eq 13.2.20: <l'm'|{n..n}_{k kappa} (chain 2,3,..k)|lm> diagonal-in-l0 check ===")
# check the l'=l+k, m'=m, kappa scalar-ish structure is consistent for k=2 vs 13.2.17 route is complex; skip.
print("  (skipped; 13.2.17/20 chain not auto-verified here)")

print()
print("=== eq 13.2.17: <l'm'|n_1mu n_1nu|lm> ===")
def me_prod_exact(lp,mp,mu,nu,l,m):
    # sum over intermediate l'' m'' of <l'm'|n_mu|l''m''><l''m''|n_nu|lm>
    tot = 0
    for lpp in (l-1,l,l+1):
        if lpp < 0: continue
        mpp = m+nu
        tot += me_exact(lp,mp,mu,lpp,mpp)*me_exact(lpp,mpp,nu,l,m)
    return simplify(tot)
def me_prod_book(lp,mp,mu,nu,l,m):
    from sympy import KroneckerDelta as KD
    kappa = mu+nu
    term1 = Rational((-1)**mu,3)*(1 if lp==l else 0)*(1 if mp==m else 0)*(1 if mu==-nu else 0)
    term2 = sqrt(Rational(2,3)*Rational(2*l+1,2*lp+1))*clebsch_gordan(l,2,lp,0,0,0)\
            *clebsch_gordan(1,1,2,mu,nu,kappa)*clebsch_gordan(l,2,lp,m,kappa,mp)
    return simplify(term1+term2)
for (lp,mp,mu,nu,l,m) in [(2,1,1,0,2,0),(3,2,1,1,1,0),(2,0,1,-1,2,1),(4,2,1,1,2,0),(1,0,1,-1,1,1)]:
    e=me_prod_exact(lp,mp,mu,nu,l,m); b=me_prod_book(lp,mp,mu,nu,l,m)
    print(f"  l'={lp} m'={mp} mu={mu:+d} nu={nu:+d} l={l} m={m}: match={simplify(e-b)==0}  (exact={e})")
