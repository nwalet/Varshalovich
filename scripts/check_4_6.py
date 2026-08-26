#!/usr/bin/env python3
r"""
Checks for Section 4.6 (Sums involving D-functions) of Chapter 4, VMK.

D^J_{MN}(a,b,g) = e^{-iMa} d^J_{MN}(b) e^{-iNg}  (wigner_d);
CG via sympy.physics.wigner.clebsch_gordan(j1,j2,j3,m1,m2,m3).
Everything is evaluated numerically (complex, ~1e-12 floor) at several angles.

  eq 4.6.1   Clebsch-Gordan series (product of two D's)
  eq 4.6.3   CG-projected sum -> {J1 J2 J} D^J
  eq 4.6.5   sum over N1,N2 -> C^{JM}_{...} D^J        (checks M_3 vs M_2)
  eq 4.6.7   double CG-projected sum -> delta {J1 J2 J}
  eq 4.6.10  product of k spin-1/2 D's
  eq 4.6.13  D^J_{MN} from Cayley-Klein a,b (eq 4.6.12/4.6.14)

Usage:  python3 check_4_6.py
"""
import os, sys, cmath, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy import I, exp, Rational, S, factorial, sqrt
from sympy.physics.wigner import clebsch_gordan as _CG
from wigner_d import wigner_d, beta as _B

fact = factorial
ANG = [(0.7, 1.1, 0.4), (1.3, 0.6, 2.0), (0.5, 2.1, 1.2)]   # (alpha,beta,gamma) floats
TOL = 1e-10

_dc = {}
def dnum(J, M, N, b):
    key = (J, M, N)
    if key not in _dc:
        _dc[key] = wigner_d(J, M, N)
    return complex(_dc[key].subs(_B, b).evalf(30))

def D(J, M, N, a, b, g):
    return cmath.exp(-1j*float(M)*a) * dnum(J, M, N, b) * cmath.exp(-1j*float(N)*g)

def CG(j1, j2, j3, m1, m2, m3):
    return complex(_CG(j1, j2, j3, m1, m2, m3))

def mrange(J):
    return [J - k for k in range(int(2*J) + 1)]

def tri(j1, j2, j3):
    return 1 if (abs(j1-j2) <= j3 <= j1+j2 and float(j1+j2+j3).is_integer()) else 0

def report(tag, worst):
    ok = worst < TOL
    print(f"  {tag:42s} {'PASS' if ok else 'FAIL'}  worst={worst:.1e}")
    return ok


def check_461(pairs):
    # D^{J1}_{M1N1} D^{J2}_{M2N2} = sum_{J,M,N} C^{JM}_{J1M1J2M2} D^J_{MN} C^{JN}_{J1N1J2N2}
    worst = 0.0
    for (J1, J2) in pairs:
        for M1 in mrange(J1):
            for N1 in mrange(J1):
                for M2 in mrange(J2):
                    for N2 in mrange(J2):
                        for (a, b, g) in ANG:
                            lhs = D(J1,M1,N1,a,b,g)*D(J2,M2,N2,a,b,g)
                            rhs = 0j
                            for J in [J1+J2-k for k in range(int(2*min(J1,J2))+1)]:
                                M, N = M1+M2, N1+N2
                                if abs(M) <= J and abs(N) <= J:
                                    rhs += CG(J1,J2,J,M1,M2,M)*D(J,M,N,a,b,g)*CG(J1,J2,J,N1,N2,N)
                            worst = max(worst, abs(lhs-rhs))
    return report("4.6.1  CG series", worst)


def check_463(triples):
    # sum_{M1M2N1N2} C^{JM}_{J1M1J2M2} D D C^{J'N}_{J1N1J2N2} = d_{JJ'}{J1J2J} D^J_{MN}
    worst = 0.0
    for (J1, J2, J) in triples:
        for Jp in [J, J+1]:
            for M in mrange(J):
                for N in mrange(J):
                    for (a, b, g) in ANG:
                        s = 0j
                        for M1 in mrange(J1):
                            M2 = M - M1
                            if abs(M2) > J2: continue
                            for N1 in mrange(J1):
                                for N2 in mrange(J2):
                                    if N1+N2 != N: continue
                                    s += (CG(J1,J2,J,M1,M2,M)*D(J1,M1,N1,a,b,g)
                                          *D(J2,M2,N2,a,b,g)*CG(J1,J2,Jp,N1,N2,N))
                        rhs = (tri(J1,J2,J)*D(J,M,N,a,b,g)) if Jp == J else 0j
                        worst = max(worst, abs(s-rhs))
    return report("4.6.3  CG-projected -> {J1J2J}D^J", worst)


def check_465(pairs):
    # sum_{N1N2} D^{J1}_{M1N1} D^{J2}_{M2N2} C^{JN}_{J1N1J2N2}
    #           = C^{JM}_{J1M1J2M2} D^J_{MN},  M=M1+M2   (M_3 must be M_2)
    worst = 0.0
    for (J1, J2) in pairs:
        for J in [J1+J2-k for k in range(int(2*min(J1,J2))+1)]:
            for M1 in mrange(J1):
                for M2 in mrange(J2):
                    M = M1+M2
                    if abs(M) > J: continue
                    for N in mrange(J):
                        for (a, b, g) in ANG:
                            s = 0j
                            for N1 in mrange(J1):
                                N2 = N - N1
                                if abs(N2) > J2: continue
                                s += D(J1,M1,N1,a,b,g)*D(J2,M2,N2,a,b,g)*CG(J1,J2,J,N1,N2,N)
                            rhs = CG(J1,J2,J,M1,M2,M)*D(J,M,N,a,b,g)
                            worst = max(worst, abs(s-rhs))
    return report("4.6.5  sum_{N1N2} -> C D^J", worst)


def check_466(pairs):
    # sum_{N1N2N} D^{J*}_{MN} D^{J1}_{M1N1} D^{J2}_{M2N2} C^{JN}_{J1N1J2N2}
    #            = C^{JM}_{J1M1J2M2}    (book sums ONLY over N1,N2,N)
    worst = 0.0
    for (J1, J2) in pairs:
        for J in [J1+J2-k for k in range(int(2*min(J1,J2))+1)]:
            for M1 in mrange(J1):
                for M2 in mrange(J2):
                    M = M1+M2
                    if abs(M) > J: continue
                    for (a, b, g) in ANG:
                        s = 0j
                        for N1 in mrange(J1):
                            for N2 in mrange(J2):
                                N = N1+N2
                                if abs(N) > J: continue
                                s += (D(J,M,N,a,b,g).conjugate()*D(J1,M1,N1,a,b,g)
                                      *D(J2,M2,N2,a,b,g)*CG(J1,J2,J,N1,N2,N))
                        worst = max(worst, abs(s - CG(J1,J2,J,M1,M2,M)))
    return report("4.6.6  sum_{N1N2N} D* D D C -> C", worst)


def check_4610(ks):
    # sum_{sum m_i=M, sum n_i=N} prod D^{1/2}_{m_i n_i} = (2J)!/sqrt(...) D^J_{MN}, J=k/2
    H = Rational(1, 2)
    worst = 0.0
    for k in ks:
        J = Rational(k, 2)
        pref0 = float(fact(2*J))
        for M in mrange(J):
            for N in mrange(J):
                for (a, b, g) in ANG:
                    s = 0j
                    for ms in itertools.product([H, -H], repeat=k):
                        if sum(ms) != M: continue
                        for ns in itertools.product([H, -H], repeat=k):
                            if sum(ns) != N: continue
                            p = 1j*0 + 1.0
                            for i in range(k):
                                p *= D(H, ms[i], ns[i], a, b, g)
                            s += p
                    denom = float(sqrt(fact(J+M)*fact(J-M)*fact(J+N)*fact(J-N)))
                    rhs = pref0/denom * D(J, M, N, a, b, g)
                    worst = max(worst, abs(s-rhs))
    return report("4.6.10 product of k spin-1/2 D's", worst)


def check_4613(Js):
    # D^J_{MN} = sqrt(...) sum_{pqrs} a^p b^q (a*)^r (-b*)^s /(p!q!r!s!),
    #   p+q+r+s=2J, p-q-r+s=2M, p+q-r-s=2N   (a=D^{1/2}_{1/2,1/2}, b=D^{1/2}_{-1/2,1/2})
    H = Rational(1, 2)
    worst = 0.0
    for (a, b, g) in ANG:
        aa = D(H, H, H, a, b, g); bb = D(H, -H, H, a, b, g)
        for J in Js:
            for M in mrange(J):
                for N in mrange(J):
                    s = 0j
                    twoJ = int(2*J)
                    for p in range(twoJ+1):
                        for q in range(twoJ+1-p):
                            for r in range(twoJ+1-p-q):
                                ss = twoJ - p - q - r
                                if ss < 0: continue
                                if p-q-r+ss != int(2*M): continue
                                if p+q-r-ss != int(2*N): continue
                                s += (aa**p * bb**q * aa.conjugate()**r * (-bb.conjugate())**ss
                                      / float(fact(p)*fact(q)*fact(r)*fact(ss)))
                    pref = float(sqrt(fact(J+M)*fact(J-M)*fact(J+N)*fact(J-N)))
                    worst = max(worst, abs(pref*s - D(J, M, N, a, b, g)))
    return report("4.6.13 Cayley-Klein a,b form of D", worst)


def main():
    print("Section 4.6 sums of D-functions\n")
    ok = True
    ok &= check_461([(Rational(1,2), S(1)), (S(1), S(1))])
    ok &= check_463([(Rational(1,2), Rational(1,2), S(1)), (S(1), S(1), S(1))])
    ok &= check_465([(Rational(1,2), S(1)), (S(1), S(1))])
    ok &= check_466([(Rational(1,2), S(1)), (S(1), S(1))])
    ok &= check_4610([2, 3, 4])
    ok &= check_4613([Rational(1,2), S(1), Rational(3,2)])
    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
