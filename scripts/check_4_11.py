#!/usr/bin/env python3
r"""
Checks for Section 4.11 (Integrals involving D-functions) of Chapter 4, VMK.

For D^J_{MN}(a,b,g)=e^{-iMa} d^J_{MN}(b) e^{-iNg}, the alpha/gamma integrals over
[0,2pi) give exact Kronecker deltas (2*pi*delta on the M-sum / N-sum); only the
beta integral  int_0^pi sin b (product of d's) db  is done numerically (quad).

  4.11.1  int D^J = d_{J0}d_{M0}d_{M'0} 8pi^2
  4.11.2  int D^{J1} D^{J2}  (product of two)
  4.11.3  int D^{J2*} D^{J1} = orthogonality
  4.11.4  int D^{J3} D^{J2} D^{J1}  (three, no conj)
  4.11.5  int D^{J3*} D^{J2} D^{J1}  (three, one conj)
  4.11.6  int sin b d^J_00 = 2 d_{J0}
  4.11.7  int sin b d d = 2/(2J+1) d_{JJ'}
  4.11.8  int sin b d^{J1} d^{J2} d^{J3} = 2/(2J3+1) C C

Usage:  python3 check_4_11.py
"""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy import Rational, S, lambdify
from sympy.physics.wigner import clebsch_gordan as _CG
from wigner_d import wigner_d, beta as _B
from scipy.integrate import quad

TOL = 1e-8
_dl = {}
def _valid(J, M): return abs(M) <= J and (J - M) == int(J - M) and (J - M) >= 0
def dfun(J, M, N):
    if not _valid(J, M) or not _valid(J, N): return None
    if (J, M, N) not in _dl:
        _dl[(J, M, N)] = lambdify(_B, wigner_d(J, M, N), "math")
    return _dl[(J, M, N)]

def betaint(fns):
    """int_0^pi sin b * prod fns(b) db (fns real-valued d-functions)."""
    def integrand(b):
        p = math.sin(b)
        for f in fns:
            p *= f(b)
        return p
    val, _ = quad(integrand, 0.0, math.pi, limit=200)
    return val

def CG(*a): return float(_CG(*a))
def mrange(J): return [J - k for k in range(int(2*J)+1)]
def kd(a, b): return 1.0 if a == b else 0.0
def report(tag, worst):
    ok = worst < TOL
    print(f"  {tag:40s} {'PASS' if ok else 'FAIL'}  worst={worst:.1e}")
    return ok

TWO_PI = 2*math.pi


def check_4111():
    # int D^J_{MM'} = delta_J0 delta_M0 delta_M'0 * 8pi^2  (J integer)
    worst = 0.0
    for J in [S(0), S(1), S(2)]:
        for M in mrange(J):
            for Mp in mrange(J):
                f = dfun(J, M, Mp)
                # alpha int: 2pi delta_{M,0}; gamma int: 2pi delta_{M',0}
                lhs = TWO_PI*kd(M,0) * TWO_PI*kd(Mp,0) * (betaint([f]) if (M==0 and Mp==0) else 0.0)
                rhs = kd(J,0)*kd(M,0)*kd(Mp,0)*8*math.pi**2
                worst = max(worst, abs(lhs-rhs))
    return report("4.11.1 int D^J", worst)


def check_4112():
    # int D^{J1}_{M1M1'} D^{J2}_{M2M2'} = (-1)^{M2-M2'} 8pi^2/(2J2+1) d_{J1J2} d_{-M1,M2} d_{-M1',M2'}
    worst = 0.0
    for J1 in [Rational(1,2), S(1)]:
        for J2 in [Rational(1,2), S(1)]:
            for M1 in mrange(J1):
                for M1p in mrange(J1):
                    for M2 in mrange(J2):
                        for M2p in mrange(J2):
                            f1, f2 = dfun(J1,M1,M1p), dfun(J2,M2,M2p)
                            lhs = TWO_PI*kd(M1+M2,0) * TWO_PI*kd(M1p+M2p,0) * (betaint([f1,f2]) if (M1+M2==0 and M1p+M2p==0) else 0.0)
                            rhs = ((-1)**int(M2-M2p) * 8*math.pi**2/(2*J2+1)
                                   * kd(J1,J2)*kd(-M1,M2)*kd(-M1p,M2p))
                            worst = max(worst, abs(lhs-rhs))
    return report("4.11.2 int D D", worst)


def check_4113():
    # int D^{J2*}_{M2M2'} D^{J1}_{M1M1'} = 8pi^2/(2J2+1) d_{J1J2} d_{M1M2} d_{M1'M2'}
    worst = 0.0
    for J1 in [Rational(1,2), S(1)]:
        for J2 in [Rational(1,2), S(1)]:
            for M1 in mrange(J1):
                for M1p in mrange(J1):
                    for M2 in mrange(J2):
                        for M2p in mrange(J2):
                            f1, f2 = dfun(J1,M1,M1p), dfun(J2,M2,M2p)
                            # D^{J2*}: e^{+iM2 a} -> alpha int 2pi delta_{M1-M2,0}
                            lhs = TWO_PI*kd(M1-M2,0)*TWO_PI*kd(M1p-M2p,0)*(betaint([f1,f2]) if (M1==M2 and M1p==M2p) else 0.0)
                            rhs = 8*math.pi**2/(2*J2+1)*kd(J1,J2)*kd(M1,M2)*kd(M1p,M2p)
                            worst = max(worst, abs(lhs-rhs))
    return report("4.11.3 int D* D (orthogonality)", worst)


def check_4114():
    # int D^{J3}_{M3M3'} D^{J2}_{M2M2'} D^{J1}_{M1M1'}
    #  = (-1)^{M3-M3'} 8pi^2/(2J3+1) C^{J3,-M3}_{J1M1J2M2} C^{J3,-M3'}_{J1M1'J2M2'}
    worst = 0.0
    tset = [(Rational(1,2),Rational(1,2),S(1)),(S(1),S(1),S(1)),(S(1),S(1),S(2))]
    for (J1,J2,J3) in tset:
        for M1 in mrange(J1):
            for M1p in mrange(J1):
                for M2 in mrange(J2):
                    for M2p in mrange(J2):
                        M3, M3p = -(M1+M2), -(M1p+M2p)   # no-conj: alpha int -> M1+M2+M3=0
                        if abs(M3) > J3 or abs(M3p) > J3: continue
                        f = [dfun(J3,M3,M3p), dfun(J2,M2,M2p), dfun(J1,M1,M1p)]
                        lhs = TWO_PI*TWO_PI*betaint(f)
                        rhs = ((-1)**int(M3-M3p)*8*math.pi**2/(2*J3+1)
                               * CG(J1,J2,J3,M1,M2,-M3)*CG(J1,J2,J3,M1p,M2p,-M3p))
                        worst = max(worst, abs(lhs-rhs))
    return report("4.11.4 int D D D", worst)


def check_4115():
    # int D^{J3*}_{M3M3'} D^{J2} D^{J1} = 8pi^2/(2J3+1) C^{J3M3}_{J1M1J2M2} C^{J3M3'}_{...}
    worst = 0.0
    tset = [(Rational(1,2),Rational(1,2),S(1)),(S(1),S(1),S(1)),(S(1),S(1),S(2))]
    for (J1,J2,J3) in tset:
        for M1 in mrange(J1):
            for M1p in mrange(J1):
                for M2 in mrange(J2):
                    for M2p in mrange(J2):
                        M3, M3p = M1+M2, M1p+M2p
                        if abs(M3) > J3 or abs(M3p) > J3: continue
                        f = [dfun(J3,M3,M3p), dfun(J2,M2,M2p), dfun(J1,M1,M1p)]
                        lhs = TWO_PI*TWO_PI*betaint(f)
                        rhs = (8*math.pi**2/(2*J3+1)
                               * CG(J1,J2,J3,M1,M2,M3)*CG(J1,J2,J3,M1p,M2p,M3p))
                        worst = max(worst, abs(lhs-rhs))
    return report("4.11.5 int D* D D", worst)


def check_411_dints():
    worst6 = 0.0
    for J in [S(0), S(1), S(2)]:
        lhs = betaint([dfun(J,0,0)]); rhs = 2*kd(J,0)
        worst6 = max(worst6, abs(lhs-rhs))
    ok = report("4.11.6 int sinb d^J_00 = 2 dJ0", worst6)
    worst7 = 0.0
    for J in [Rational(1,2),S(1),Rational(3,2)]:
        for Jp in [J, J+1]:
            for M in mrange(J):
                for Mp in mrange(J):
                    fj = dfun(J,M,Mp); fjp = dfun(Jp,M,Mp)
                    lhs = betaint([fj,fjp]) if fjp else 0.0
                    rhs = 2/(2*J+1)*kd(J,Jp)
                    worst7 = max(worst7, abs(lhs-rhs))
    ok &= report("4.11.7 int sinb d d = 2/(2J+1) dJJ'", worst7)
    worst8 = 0.0
    for (J1,J2,J3) in [(Rational(1,2),Rational(1,2),S(1)),(S(1),S(1),S(2))]:
        for M1 in mrange(J1):
            for M1p in mrange(J1):
                for M2 in mrange(J2):
                    for M2p in mrange(J2):
                        M3, M3p = M1+M2, M1p+M2p
                        if abs(M3)>J3 or abs(M3p)>J3: continue
                        lhs = betaint([dfun(J1,M1,M1p),dfun(J2,M2,M2p),dfun(J3,M3,M3p)])
                        rhs = 2/(2*J3+1)*CG(J1,J2,J3,M1,M2,M3)*CG(J1,J2,J3,M1p,M2p,M3p)
                        worst8 = max(worst8, abs(lhs-rhs))
    ok &= report("4.11.8 int sinb d d d = 2/(2J3+1) C C", worst8)
    return ok


def main():
    print("Section 4.11 integrals of D-functions\n")
    ok = True
    ok &= check_4111(); ok &= check_4112(); ok &= check_4113()
    ok &= check_4114(); ok &= check_4115(); ok &= check_411_dints()
    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
