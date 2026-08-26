#!/usr/bin/env python3
r"""
Checks for Section 4.5 (Rotation matrix U^J_{MM'}(omega; Theta, Phi)) of
Chapter 4, Varshalovich-Moskalev-Khersonskii.

Trusted reference: u_function() in Ufunction.py -- the Cayley-Klein closed
form (eq 4.5.4), validated there against the matrix exponential
exp(-i omega n.J).  D^J_{MM'}(a,b,g) = e^{-iMa} d^J_{MM'}(b) e^{-iM'g}
(wigner_d).  Everything is compared numerically at several (omega,Theta,Phi).

  eq 4.5.3   U as a product of two D-functions
  eq 4.5.4   Cayley-Klein sum form (both M+M' cases), as printed
  eq 4.5.6   U via d^J_{MM'}(xi),  sin(xi/2)=sin(w/2)sin(Theta)
  eq 4.5.17  inverse:  U(-w;Th,Ph) = U(w;pi-Th,pi+Ph)
  eq 4.5.18  conjugation
  eq 4.5.19  argument-sign reversal (3)
  eq 4.5.20  periodicity (3)
  eq 4.5.21  half-period shifts (6)
  eq 4.5.22  M<->M' permutation (2)
  eq 4.5.28-30, 4.5.32  special cases

Usage:  python3 check_4_5.py
"""
import os, sys, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy import (I, exp, pi, Rational, S, Abs, sqrt, factorial, cos, sin,
                   tan, asin, conjugate)
from Ufunction import u_function, omega as _W, Theta as _TH, Phi as _PH
from wigner_d import wigner_d, beta as _B

fact = factorial
# sample (omega, Theta, Phi) with Theta,Phi away from 0/pi (u_function is
# singular at Theta=0) and generic to avoid coincidences
PTS = [(Rational(1,3)*pi, Rational(1,4)*pi, Rational(1,5)*pi),
       (Rational(2,5)*pi, Rational(1,3)*pi, Rational(2,7)*pi),
       (Rational(1,4)*pi, Rational(2,5)*pi, Rational(1,6)*pi)]
TOL = 1e-11   # values via complex() -> double precision

_uc = {}
def Uref(J, M, Mp, w, th, ph):
    key = (J, M, Mp)
    if key not in _uc:
        _uc[key] = u_function(J, M, Mp)          # expr in _W,_TH,_PH
    return complex(_uc[key].subs({_W: w, _TH: th, _PH: ph}).evalf(30))

_dc = {}
def dref(J, M, Mp, b):
    key = (J, M, Mp)
    if key not in _dc:
        _dc[key] = wigner_d(J, M, Mp)
    return _dc[key].subs(_B, b)

def Dv(J, M, Mp, a, b, g):
    return complex((exp(-I*M*a)*dref(J, M, Mp, b)*exp(-I*Mp*g)).evalf(30))

def jrange(J):
    return [J - k for k in range(int(2*J) + 1)]

def ph_(x):
    return complex(exp(I*pi*x).evalf(30))


def run(tag, lhs, rhs, Jset, pts=PTS):
    worst = 0.0
    for J in Jset:
        for M in jrange(J):
            for Mp in jrange(J):
                for (w, th, pp) in pts:
                    worst = max(worst, abs(lhs(J,M,Mp,w,th,pp) - rhs(J,M,Mp,w,th,pp)))
    ok = worst < TOL
    print(f"  {tag:56s} {'PASS' if ok else 'FAIL'}  worst={worst:.1e}")
    return ok


# eq 4.5.3  U = sum_{M''} D(Ph,Th,-Ph) e^{-iM'' w} D(Ph,-Th,-Ph)
def U_453(J, M, Mp, w, th, pp):
    return sum(Dv(J,M,Mpp,pp,th,-pp) * cmath.exp(-1j*float(Mpp)*float(w))
               * Dv(J,Mpp,Mp,pp,-th,-pp) for Mpp in jrange(J))

# eq 4.5.4  as printed in Chap4.tex
def U_454(J, M, Mp, w, th, pp):
    v = sin(w/2)*sin(th); u = cos(w/2) - I*sin(w/2)*cos(th)
    pre = sqrt(fact(J+M)*fact(J-M)*fact(J+Mp)*fact(J-Mp))
    if M + Mp >= 0:
        base = (-I*v)**(2*J) * (u/(-I*v))**(M+Mp) * exp(-I*(M-Mp)*pp)
        ss = [s for s in range(0, 40)
              if min(s, s+M+Mp, J-M-s, J-Mp-s) >= 0]
        tot = sum((1 - v**-2)**s / (fact(s)*fact(s+M+Mp)*fact(J-M-s)*fact(J-Mp-s)) for s in ss)
    else:
        u2 = cos(w/2) + I*sin(w/2)*cos(th)          # u*
        base = (-I*v)**(2*J) * (u2/(-I*v))**(-(M+Mp)) * exp(-I*(M-Mp)*pp)
        ss = [s for s in range(0, 40)
              if min(s, s-M-Mp, J+M-s, J+Mp-s) >= 0]
        tot = sum((1 - v**-2)**s / (fact(s)*fact(s-M-Mp)*fact(J+M-s)*fact(J+Mp-s)) for s in ss)
    return complex((pre*base*tot).evalf(30))

# eq 4.5.6  U = i^{M-M'} e^{-i(M-M')Ph} ((1-i tan(w/2)cosTh)/sqrt(1+tan^2 cos^2))^{M+M'} d(xi)
def U_456(J, M, Mp, w, th, pp):
    xi = 2*asin(sin(w/2)*sin(th))
    fac = (1 - I*tan(w/2)*cos(th)) / sqrt(1 + tan(w/2)**2*cos(th)**2)
    val = (I**(M-Mp) * exp(-I*(M-Mp)*pp) * fac**(M+Mp) * dref(J,M,Mp,xi))
    return complex(val.evalf(30))


def main():
    Jset = [Rational(1,2), S(1), Rational(3,2)]
    print(f"Section 4.5 U-matrix  (J up to {Jset[-1]})\n")
    ok = True
    print("explicit forms")
    ok &= run("4.5.3  D-product form", U_453, Uref, Jset)
    ok &= run("4.5.4  Cayley-Klein sum (as printed)", U_454, Uref, Jset)
    ok &= run("4.5.6  d(xi) form", U_456, Uref, Jset)

    print("\nproperties")
    ok &= run("4.5.17 U(-w;Th,Ph)=U(w;pi-Th,pi+Ph)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,-w,th,pp),
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,pi-th,pi+pp), Jset)
    ok &= run("4.5.18a U*=U_{M'M}(w;pi-Th,pi+Ph)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,th,pp).conjugate(),
              lambda J,M,Mp,w,th,pp: Uref(J,Mp,M,w,pi-th,pi+pp), Jset)
    ok &= run("4.5.18b U*=U_{M'M}(-w;Th,Ph)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,th,pp).conjugate(),
              lambda J,M,Mp,w,th,pp: Uref(J,Mp,M,-w,th,pp), Jset)
    ok &= run("4.5.19a U(-w;Th,Ph)=(-1)^{M-M'}U_{-M',-M}(w;Th,Ph)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,-w,th,pp),
              lambda J,M,Mp,w,th,pp: ph_(M-Mp)*Uref(J,-Mp,-M,w,th,pp), Jset)
    ok &= run("4.5.19b U(w;-Th,Ph)=(-1)^{M-M'}U(w;Th,Ph)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,-th,pp),
              lambda J,M,Mp,w,th,pp: ph_(M-Mp)*Uref(J,M,Mp,w,th,pp), Jset)
    ok &= run("4.5.19c U(w;Th,-Ph)=U_{M'M}(w;Th,Ph)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,th,-pp),
              lambda J,M,Mp,w,th,pp: Uref(J,Mp,M,w,th,pp), Jset)
    ok &= run("4.5.20a U(w+4pi)=U(w)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w+4*pi,th,pp),
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,th,pp), Jset)
    ok &= run("4.5.21a U(w+2pi)=(-1)^{2J}U(w)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w+2*pi,th,pp),
              lambda J,M,Mp,w,th,pp: ph_(2*J)*Uref(J,M,Mp,w,th,pp), Jset)
    ok &= run("4.5.21b U(2pi-w)=(-1)^{M+M'}U_{-M',-M}(w)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,2*pi-w,th,pp),
              lambda J,M,Mp,w,th,pp: ph_(M+Mp)*Uref(J,-Mp,-M,w,th,pp), Jset)
    ok &= run("4.5.21c U(w;Th+pi,Ph)=(-1)^{M-M'}U_{-M',-M}(w;Th,Ph)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,th+pi,pp),
              lambda J,M,Mp,w,th,pp: ph_(M-Mp)*Uref(J,-Mp,-M,w,th,pp), Jset)
    ok &= run("4.5.21d U(w;pi-Th,Ph)=U_{-M',-M}(w;Th,Ph)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,pi-th,pp),
              lambda J,M,Mp,w,th,pp: Uref(J,-Mp,-M,w,th,pp), Jset)
    ok &= run("4.5.21e U(w;Th,Ph+pi)=(-1)^{M-M'}U(w;Th,Ph)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,th,pp+pi),
              lambda J,M,Mp,w,th,pp: ph_(M-Mp)*Uref(J,M,Mp,w,th,pp), Jset)
    ok &= run("4.5.21f U(w;Th,pi-Ph)=(-1)^{M-M'}U_{M'M}(w;Th,Ph)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,th,pi-pp),
              lambda J,M,Mp,w,th,pp: ph_(M-Mp)*Uref(J,Mp,M,w,th,pp), Jset)
    ok &= run("4.5.22a U_{M'M}=U_{MM'}(w;Th,-Ph)",
              lambda J,M,Mp,w,th,pp: Uref(J,Mp,M,w,th,pp),
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,th,-pp), Jset)
    ok &= run("4.5.22b U_{-M-M'}=(-1)^{M-M'}U_{MM'}(w;pi-Th,pi-Ph)",
              lambda J,M,Mp,w,th,pp: Uref(J,-M,-Mp,w,th,pp),
              lambda J,M,Mp,w,th,pp: ph_(M-Mp)*Uref(J,M,Mp,w,pi-th,pi-pp), Jset)

    print("\nspecial cases")
    ok &= run("4.5.28 U(w;pi/2,0)=(i)^{M-M'}d(w) [corrected]",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,pi/2,S(0)),
              lambda J,M,Mp,w,th,pp: complex(((I)**(M-Mp)*dref(J,M,Mp,w)).evalf(30)), Jset)
    ok &= run("4.5.29 U(w;pi/2,pi/2)=d(w)",
              lambda J,M,Mp,w,th,pp: Uref(J,M,Mp,w,pi/2,pi/2),
              lambda J,M,Mp,w,th,pp: complex(dref(J,M,Mp,w).evalf(30)), Jset)
    ok &= run("4.5.30 U(w;0,Ph)=delta e^{-iMw}  [via 4.5.3]",
              lambda J,M,Mp,w,th,pp: U_453(J,M,Mp,w,S(0),pp),
              lambda J,M,Mp,w,th,pp: (cmath.exp(-1j*float(M)*float(w)) if M==Mp else 0j), Jset)
    # 4.5.32 corner elements
    ok &= run("4.5.32 U_{JJ}=(cos w/2 - i sin w/2 cosTh)^{2J}",
              lambda J,M,Mp,w,th,pp: Uref(J,J,J,w,th,pp),
              lambda J,M,Mp,w,th,pp: complex(((cos(w/2)-I*sin(w/2)*cos(th))**(2*J)).evalf(30)), Jset)
    ok &= run("4.5.32 U_{J,-J}=(-i sin w/2 sinTh e^{-iPh})^{2J}",
              lambda J,M,Mp,w,th,pp: Uref(J,J,-J,w,th,pp),
              lambda J,M,Mp,w,th,pp: complex(((-I*sin(w/2)*sin(th)*exp(-I*pp))**(2*J)).evalf(30)), Jset)

    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
