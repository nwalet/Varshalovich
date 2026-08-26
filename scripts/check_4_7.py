#!/usr/bin/env python3
r"""
Checks for Section 4.7 (Addition of rotations) of Chapter 4, VMK.

d^J from wigner_d;  D^J_{MN}(a,b,g)=e^{-iMa}d^J_{MN}(b)e^{-iNg};  CG from sympy.
Numeric (complex) checks at sample angles.

  eq 4.7.4    unitarity: sum_{M''} D_{MM''} D*_{M'M''} = delta_{MM'}
  eq 4.7.7    d-addition, phi=0, b1+b2<=pi:  sum d(b1)d(b2)=d(b1+b2)
  eq 4.7.8    d-addition, phi=0, b1+b2>pi:   ... =(-1)^{M+M'} d(2pi-b1-b2)
  eq 4.7.9    d-addition, phi=pi, b1>=b2:     sum (-1)^{M''-M'} d d = d(b1-b2)
  eq 4.7.10   b1=b2 case of 4.7.9 -> delta_{MM'}
  eq 4.7.15   identical rotations, a=g=0 -> d^j_{mm'}(2b)
  eq 4.7.17   character sum chi^J(R1R2)=sin[(2J+1)w/2]/sin(w/2), w from 4.7.18
  eq 4.7.19   character sum chi^J(R1 R2^{-1}), w' from 4.7.20

Usage:  python3 check_4_7.py
"""
import os, sys, cmath, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy import Rational, S
from sympy.physics.wigner import clebsch_gordan as _CG
from wigner_d import wigner_d, beta as _B

TOL = 1e-10
_dc = {}
def dnum(J, M, N, b):
    key = (J, M, N)
    if key not in _dc:
        _dc[key] = wigner_d(J, M, N)
    return complex(_dc[key].subs(_B, float(b)).evalf(30))

def D(J, M, N, a, b, g):
    return cmath.exp(-1j*float(M)*a)*dnum(J, M, N, b)*cmath.exp(-1j*float(N)*g)

def CG(j1, j2, j3, m1, m2, m3):
    return complex(_CG(j1, j2, j3, m1, m2, m3))

def mrange(J):
    return [J - k for k in range(int(2*J) + 1)]

def report(tag, worst):
    ok = worst < TOL
    print(f"  {tag:46s} {'PASS' if ok else 'FAIL'}  worst={worst:.1e}")
    return ok


def check_474(Js):
    worst = 0.0
    for J in Js:
        for M in mrange(J):
            for Mp in mrange(J):
                for (a, b, g) in [(0.7,1.1,0.4),(1.3,0.6,2.0)]:
                    s1 = sum(D(J,M,Mpp,a,b,g)*D(J,Mp,Mpp,a,b,g).conjugate() for Mpp in mrange(J))
                    s2 = sum(D(J,M,Mpp,a,b,g)*D(J,Mpp,Mp,-g,-b,-a) for Mpp in mrange(J))
                    tgt = 1.0 if M == Mp else 0.0
                    worst = max(worst, abs(s1-tgt), abs(s2-tgt))
    return report("4.7.4  unitarity", worst)


def _dadd(Js, b1, b2, phase, rhs_fn):
    worst = 0.0
    for J in Js:
        for M in mrange(J):
            for Mp in mrange(J):
                s = sum(phase(Mpp, Mp)*dnum(J,M,Mpp,b1)*dnum(J,Mpp,Mp,b2) for Mpp in mrange(J))
                worst = max(worst, abs(s - rhs_fn(J, M, Mp)))
    return worst

def check_477(Js):
    b1, b2 = 0.5, 0.7      # b1+b2 < pi
    w = _dadd(Js, b1, b2, lambda Mpp, Mp: 1.0, lambda J,M,Mp: dnum(J,M,Mp,b1+b2))
    return report("4.7.7  phi=0, b1+b2<=pi", w)

def check_478(Js):
    b1, b2 = 1.8, 1.9      # b1+b2 > pi
    w = _dadd(Js, b1, b2, lambda Mpp, Mp: 1.0,
              lambda J,M,Mp: complex((-1)**int(M+Mp))*dnum(J,M,Mp,2*math.pi-b1-b2))
    return report("4.7.8  phi=0, b1+b2>pi", w)

def check_479(Js):
    b1, b2 = 1.5, 0.6      # b1 >= b2
    w = _dadd(Js, b1, b2, lambda Mpp, Mp: complex((-1)**int(Mpp-Mp)),
              lambda J,M,Mp: dnum(J,M,Mp,b1-b2))
    return report("4.7.9  phi=pi, b1>=b2", w)

def check_4710(Js):
    b = 0.9
    w = _dadd(Js, b, b, lambda Mpp, Mp: complex((-1)**int(Mpp-Mp)),
              lambda J,M,Mp: 1.0 if M == Mp else 0.0)
    return report("4.7.10 b1=b2 -> delta", w)


def check_4715(js):
    # a=g=0: sum_J sum_{m''} C^{J,m+m''}_{jm,jm''} d^J_{m+m'',m''+m'}(b) C^{J,m''+m'}_{jm'',jm'}
    #        = d^j_{mm'}(2b)
    worst = 0.0
    b = 0.8
    for j in js:
        for m in mrange(j):
            for mp in mrange(j):
                s = 0j
                for mpp in mrange(j):
                    Mtop, Nbot = m+mpp, mpp+mp
                    for J in [2*j - k for k in range(int(2*j)+1)]:
                        if abs(Mtop) > J or abs(Nbot) > J: continue
                        s += (CG(j,j,J,m,mpp,Mtop)*dnum(J,Mtop,Nbot,b)
                              *CG(j,j,J,mpp,mp,Nbot))
                worst = max(worst, abs(s - dnum(j,m,mp,2*b)))
    return report("4.7.15 identical rotations (a=g=0)", worst)


# ---- character sums 4.7.17-4.7.20 ----
def chi_prod(J, R1, R2, conj2):
    a1,b1,g1 = R1; a2,b2,g2 = R2
    s = 0j
    for M in mrange(J):
        for Mp in mrange(J):
            d2 = D(J,M,Mp,a2,b2,g2) if conj2 else D(J,Mp,M,a2,b2,g2)
            if conj2:
                d2 = d2.conjugate()
            s += D(J,M,Mp,a1,b1,g1)*d2
    return s

def chi_formula(J, half):   # sin[(2J+1)*half]/sin(half), half = w/2
    num = math.sin((2*float(J)+1)*half)
    return num/math.sin(half)

def check_4717(Js):
    R1 = (0.6, 1.0, 0.5); R2 = (1.2, 0.7, 0.9)
    a1,b1,g1 = R1; a2,b2,g2 = R2
    # eq 4.7.18 two expressions for cos(w/2)
    c1 = (math.cos(b1/2)*math.cos(b2/2)*math.cos((a1+g1+a2+g2)/2)
          - math.sin(b1/2)*math.sin(b2/2)*math.cos((a1-g1-a2+g2)/2))
    c2 = (math.cos((b1+b2)/2)*math.cos((a1+g2)/2)*math.cos((g1+a2)/2)
          - math.cos((b1-b2)/2)*math.sin((a1+g2)/2)*math.sin((g1+a2)/2))
    worst_c = abs(c1-c2)
    half = math.acos(max(-1.0, min(1.0, c1)))
    worst = worst_c
    for J in Js:
        lhs = chi_prod(J, R1, R2, conj2=False)   # sum D(R1)_{MM'} D(R2)_{M'M}
        worst = max(worst, abs(lhs.imag), abs(lhs.real - chi_formula(J, half)))
    print(f"    (4.7.18 cos(w/2) two forms agree: {worst_c:.1e})")
    return report("4.7.17 chi(R1R2) + 4.7.18 angle", worst)

def check_4719(Js):
    R1 = (0.6, 1.0, 0.5); R2 = (1.2, 0.7, 0.9)
    a1,b1,g1 = R1; a2,b2,g2 = R2
    c1 = (math.cos(b1/2)*math.cos(b2/2)*math.cos((a1+g1-a2-g2)/2)
          + math.sin(b1/2)*math.sin(b2/2)*math.cos((a1-g1-a2+g2)/2))
    c2 = (math.cos((b1-b2)/2)*math.cos((a1-a2)/2)*math.cos((g1-g2)/2)
          - math.cos((b1+b2)/2)*math.sin((a1-a2)/2)*math.sin((g1-g2)/2))
    worst_c = abs(c1-c2)
    half = math.acos(max(-1.0, min(1.0, c1)))
    worst = worst_c
    for J in Js:
        lhs = chi_prod(J, R1, R2, conj2=True)     # sum D(R1)_{MM'} D*(R2)_{MM'}
        worst = max(worst, abs(lhs.imag), abs(lhs.real - chi_formula(J, half)))
    print(f"    (4.7.20 cos(w'/2) two forms agree: {worst_c:.1e})")
    return report("4.7.19 chi(R1 R2^-1) + 4.7.20 angle", worst)


def main():
    print("Section 4.7 addition of rotations\n")
    Js = [Rational(1,2), S(1), Rational(3,2), S(2)]
    js = [Rational(1,2), S(1), Rational(3,2)]
    ok = True
    ok &= check_474(Js)
    ok &= check_477(Js); ok &= check_478(Js); ok &= check_479(Js); ok &= check_4710(Js)
    ok &= check_4715(js)
    ok &= check_4717(Js); ok &= check_4719(Js)
    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
