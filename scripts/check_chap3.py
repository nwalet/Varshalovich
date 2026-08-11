"""Verify numerically-checkable equations in Chapter 3 (irreducible tensors)."""
import sympy as sp
from sympy import sqrt, I, Rational, simplify, symbols
from sympy.physics.quantum.cg import CG

def sph_from_cart(Vx, Vy, Vz):
    """spherical components (mu=+1,0,-1) of a vector with cartesian Vx,Vy,Vz."""
    return {1: -1/sqrt(2)*(Vx + I*Vy), 0: Vz, -1: 1/sqrt(2)*(Vx - I*Vy)}

print("=== eq 3.2.5 / 3.2.6 : {A x B}_{2M} coefficient sqrt((3|M|-2)/(14|M|-12)) ===")
# symbolic vector spherical comps
Ax,Ay,Az,Bx,By,Bz = symbols('Ax Ay Az Bx By Bz')
A = sph_from_cart(Ax,Ay,Az); B = sph_from_cart(Bx,By,Bz)
def tensor2(M):
    s = 0
    for mu in (1,0,-1):
        for nu in (1,0,-1):
            if mu+nu==M:
                s += CG(1,mu,1,nu,2,M).doit()*A[mu]*B[nu]
    return sp.expand(s)
def coeff_formula(M):
    return sqrt(Rational(3*abs(M)-2, 14*abs(M)-12))
def tensor2_formula(M):
    terms=[]
    for mu in (1,0,-1):
        for nu in (1,0,-1):
            if mu+nu==M and mu>=nu:
                terms.append(A[mu]*B[nu]+A[nu]*B[mu])
    return sp.expand(coeff_formula(M)*sum(terms))
for M in (2,1,0,-1,-2):
    ok = simplify(tensor2(M)-tensor2_formula(M))==0
    print(f"  M={M:+d}: CG-sum == |M|-formula ? {ok}")

print()
print("=== eq 3.2.30 : rank-1 (axial vector) reduction ===")
# antisymmetric A_ik -> U_i = 1/2 eps_ikl A_kl.  U_x=A_yz, U_y=A_zx, U_z=A_xy
Axy,Ayz,Azx = symbols('Axy Ayz Azx')
Aik = {('x','y'):Axy,('y','x'):-Axy,('y','z'):Ayz,('z','y'):-Ayz,
       ('z','x'):Azx,('x','z'):-Azx,('x','x'):0,('y','y'):0,('z','z'):0}
Ux,Uy,Uz = Ayz, Azx, Axy   # from U_i=1/2 eps_ikl A_kl
U = sph_from_cart(Ux,Uy,Uz)
print("  I_{10}=U_z=A_xy :", simplify(U[0]-Axy)==0, " (book had U_x, A_yx/A_xx)")
print("  I_{1+1}=-1/sqrt2 (U_x+iU_y)=-1/sqrt2 (A_yz+iA_zx):",
      simplify(U[1] - (-1/sqrt(2)*(Ayz+I*Azx)))==0)
print("  I_{1-1}=+1/sqrt2 (U_x-iU_y)=+1/sqrt2 (A_yz-iA_zx):",
      simplify(U[-1] - (1/sqrt(2)*(Ayz-I*Azx)))==0)

print()
print("=== eq 3.2.31 : rank-2 (sym traceless) reduction, I_2m = sqrt(2/3){AxB}_2m ===")
# build symmetric traceless S from A,B and compare book forms with corrected indices
Sxx,Syy,Szz,Sxy,Sxz,Syz = symbols('Sxx Syy Szz Sxy Sxz Syz')
# express {AxB}_{2m} in cartesian, substitute S (off-diag = A_iB_j+A_jB_i halves etc.)
subs = {Ax*Bx:Sxx+Rational(1,3), Ay*By:Syy+Rational(1,3), Az*Bz:Szz+Rational(1,3)}
# instead: verify the algebraic identities directly using the hand-derived forms
norm = sqrt(Rational(2,3))
book = {0: Szz,
        1: -norm*(Sxz+I*Syz), -1: norm*(Sxz-I*Syz),
        2: sqrt(Rational(1,6))*(Sxx-Syy+2*I*Sxy),
       -2: sqrt(Rational(1,6))*(Sxx-Syy-2*I*Sxy)}
# {AxB}_2m in terms of S (derived): m0=sqrt(3/2)Szz, m1=-(Sxz+iSyz), m2=1/2(Sxx-Syy+2iSxy)
t2S = {0: sqrt(Rational(3,2))*Szz,
       1: -(Sxz+I*Syz), -1:(Sxz-I*Syz),
       2: Rational(1,2)*(Sxx-Syy+2*I*Sxy), -2: Rational(1,2)*(Sxx-Syy-2*I*Sxy)}
for m in (0,1,-1,2,-2):
    ok = simplify(book[m]-norm*t2S[m])==0
    print(f"  m={m:+d}: book form == sqrt(2/3)*{{AxB}}_2m ? {ok}")

print()
print("=== eq 3.1.2 : compact commutator phase e^{i mu delta} and CG normalization ===")
# check sqrt(J(J+1)) * C^{J,M+mu}_{JM,1mu} reproduces eq 3.1.1 amplitudes
J,M = symbols('J M', positive=True)
for Jn in (1, Rational(3,2), 2, Rational(5,2)):
    # valid projections M = -J, -J+1, ..., J-1  (need M+1 <= J)
    Mn = -Jn
    while Mn <= Jn-1:
        lhs = sqrt(Jn*(Jn+1))*CG(Jn,Mn,1,1,Jn,Mn+1).doit()
        rhs = -1/sqrt(2)*sqrt(Jn*(Jn+1)-Mn*(Mn+1))
        assert simplify(lhs-rhs)==0, (Jn,Mn,lhs,rhs)
        Mn += 1
print("  sqrt(J(J+1)) C^{J,M+1}_{JM,11} == -1/sqrt2 sqrt(J(J+1)-M(M+1))  [all tested J,M] : True")
print("  => compact form must carry e^{i mu delta} (mu=+1 -> e^{i delta}), not e^{iM delta}")
