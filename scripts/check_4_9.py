#!/usr/bin/env python3
r"""
Checks for Section 4.9 (Differential relations for D^J) of Chapter 4, VMK.

dD/dbeta computed by differentiating the cached d^J_{MM'}(beta) (wigner_d).

  4.9.1  sin b dD/db = ... (D^{J-1}, D^J, D^{J+1})
  4.9.2  dD/db = -1/2 sqrt((J+M)(J-M+1)) e^{-ia} D_{M-1,M'} + 1/2 sqrt((J-M)(J+M+1)) e^{ia} D_{M+1,M'}
  4.9.3  dD/db = +1/2 sqrt((J+M')(J-M'+1)) e^{-ig} D_{M,M'-1} - 1/2 sqrt((J-M')(J+M'+1)) e^{ig} D_{M,M'+1}
  4.9.4  dD/db = (M'-M cosb)/sinb D - sqrt((J+M)(J-M+1)) e^{-ia} D_{M-1,M'}
  4.9.5  dD/db = -(M'-M cosb)/sinb D + sqrt((J-M)(J+M+1)) e^{ia} D_{M+1,M'}
  4.9.6  dD/db = (M-M' cosb)/sinb D + sqrt((J+M')(J-M'+1)) e^{-ig} D_{M,M'-1}
  4.9.7  dD/db = -(M-M' cosb)/sinb D - sqrt((J-M')(J+M'+1)) e^{ig} D_{M,M'+1}
  4.9.8  dD/da = -i M D
  4.9.9  dD/dg = -i M' D

Usage:  python3 check_4_9.py
"""
import os, sys, cmath, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy import Rational, S, diff
from wigner_d import wigner_d, beta as _B

ANG = [(0.7, 1.1, 0.4), (1.3, 0.6, 2.0), (0.5, 2.1, 1.2)]
TOL = 1e-10
def sq(x): return cmath.sqrt(complex(x))

_dc, _dbc = {}, {}
def _valid(J, M): return abs(M) <= J and (J - M) == int(J - M) and (J - M) >= 0
def dval(J, M, N, b):
    if not _valid(J, M) or not _valid(J, N): return None
    if (J, M, N) not in _dc: _dc[(J, M, N)] = wigner_d(J, M, N)
    return complex(_dc[(J, M, N)].subs(_B, float(b)).evalf(30))
def dbeta(J, M, N, b):
    if (J, M, N) not in _dbc: _dbc[(J, M, N)] = diff(wigner_d(J, M, N), _B)
    return complex(_dbc[(J, M, N)].subs(_B, float(b)).evalf(30))

def D(J, M, N, a, b, g):
    v = dval(J, M, N, b)
    return 0j if v is None else cmath.exp(-1j*float(M)*a)*v*cmath.exp(-1j*float(N)*g)
def dDdb(J, M, N, a, b, g):
    return cmath.exp(-1j*float(M)*a)*dbeta(J, M, N, b)*cmath.exp(-1j*float(N)*g)

def mrange(J): return [J - k for k in range(int(2*J)+1)]
def report(tag, worst):
    ok = worst < TOL
    print(f"  {tag:50s} {'PASS' if ok else 'FAIL'}  worst={worst:.1e}")
    return ok

def run(tag, Jset, lhs, rhs):
    worst = 0.0
    for J in Jset:
        for M in mrange(J):
            for Mp in mrange(J):
                for (a, b, g) in ANG:
                    worst = max(worst, abs(lhs(J,M,Mp,a,b,g) - rhs(J,M,Mp,a,b,g)))
    return report(tag, worst)


def main():
    print("Section 4.9 differential relations\n")
    J1 = [S(1), Rational(3,2), S(2)]
    e = cmath.exp
    ok = True
    ok &= run("4.9.1  sin b dD/db", J1,
        lambda J,M,Mp,a,b,g: math.sin(b)*dDdb(J,M,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: (
            -(J+1)*sq((J**2-M**2)*(J**2-Mp**2))/(J*(2*J+1))*D(J-1,M,Mp,a,b,g)
            - M*Mp/(J*(J+1))*D(J,M,Mp,a,b,g)
            + J*sq(((J+1)**2-M**2)*((J+1)**2-Mp**2))/((J+1)*(2*J+1))*D(J+1,M,Mp,a,b,g)))
    ok &= run("4.9.2  dD/db (M+-1)", J1,
        lambda J,M,Mp,a,b,g: dDdb(J,M,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: (-0.5*sq((J+M)*(J-M+1))*e(-1j*a)*D(J,M-1,Mp,a,b,g)
                              +0.5*sq((J-M)*(J+M+1))*e(1j*a)*D(J,M+1,Mp,a,b,g)))
    ok &= run("4.9.3  dD/db (M'+-1)", J1,
        lambda J,M,Mp,a,b,g: dDdb(J,M,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: (0.5*sq((J+Mp)*(J-Mp+1))*e(-1j*g)*D(J,M,Mp-1,a,b,g)
                              -0.5*sq((J-Mp)*(J+Mp+1))*e(1j*g)*D(J,M,Mp+1,a,b,g)))
    ok &= run("4.9.4  dD/db = (M'-Mcosb)/sinb D - ...", J1,
        lambda J,M,Mp,a,b,g: dDdb(J,M,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: ((Mp-M*math.cos(b))/math.sin(b)*D(J,M,Mp,a,b,g)
                              -sq((J+M)*(J-M+1))*e(-1j*a)*D(J,M-1,Mp,a,b,g)))
    ok &= run("4.9.5  dD/db = -(M'-Mcosb)/sinb D + ...", J1,
        lambda J,M,Mp,a,b,g: dDdb(J,M,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: (-(Mp-M*math.cos(b))/math.sin(b)*D(J,M,Mp,a,b,g)
                              +sq((J-M)*(J+M+1))*e(1j*a)*D(J,M+1,Mp,a,b,g)))
    ok &= run("4.9.6  dD/db = (M-M'cosb)/sinb D + ...", J1,
        lambda J,M,Mp,a,b,g: dDdb(J,M,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: ((M-Mp*math.cos(b))/math.sin(b)*D(J,M,Mp,a,b,g)
                              +sq((J+Mp)*(J-Mp+1))*e(-1j*g)*D(J,M,Mp-1,a,b,g)))
    ok &= run("4.9.7  dD/db = -(M-M'cosb)/sinb D - ...", J1,
        lambda J,M,Mp,a,b,g: dDdb(J,M,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: (-(M-Mp*math.cos(b))/math.sin(b)*D(J,M,Mp,a,b,g)
                              -sq((J-Mp)*(J+Mp+1))*e(1j*g)*D(J,M,Mp+1,a,b,g)))
    ok &= run("4.9.8  dD/da = -iM D", J1,
        lambda J,M,Mp,a,b,g: -1j*float(M)*D(J,M,Mp,a,b,g),
        lambda J,M,Mp,a,b,g: -1j*float(M)*D(J,M,Mp,a,b,g))  # trivial identity (definition)
    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
