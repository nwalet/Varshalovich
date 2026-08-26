#!/usr/bin/env python3
r"""
Checks for Section 4.3 (Explicit forms of the Wigner D-functions) of Chapter 4,
Varshalovich-Moskalev-Khersonskii.

Reference d^J_{M M'}(beta) is the validated VMK-convention helper in
``wigner_d.py`` (itself proof-by-points against the Wigner sum formula).
Every form below is compared to it numerically at several interior angles
beta in (0, pi), to 30-digit precision, over a range of (J, M, M').

  Sec 4.3.1  trig sums
    eq 4.3.2   (-1)^{J-M'} form, cos^{M+M'+2k}
    eq 4.3.3   (-1)^{J+M'} form, cos^{2k-M-M'}
    eq 4.3.4   plain form,       cos^{2J-2k+M-M'}
    eq 4.3.5   (-1)^{M-M'} form, cos^{2J-2k-M+M'}
    eq 4.3.6   general Clebsch-Gordan form (J1, J2 free with J1-J2=M')

  Sec 4.3.2  differential (Rodrigues) representations
    eq 4.3.7 .. 4.3.10
  Sec 4.3.4  Jacobi-polynomial form
    eq 4.3.13 (with mu,nu,s,xi from 4.3.14/4.3.15)
  Sec 4.3.5  hypergeometric forms
    eq 4.3.16 .. 4.3.23

Usage:  python3 check_4_3.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sympy import (Rational, S, sqrt, factorial, cos, sin, symbols, pi,
                   diff, simplify, Abs, jacobi, nsimplify, Integer, prod)
from wigner_d import wigner_d


def _poch(a, k):
    return prod([a + i for i in range(k)]) if k > 0 else S(1)


def hyper(top, bot, z):
    """Terminating generalized hypergeometric (some top param is a
    non-positive integer, so the series is a finite polynomial in z).
    Signature mirrors sympy.hyper: hyper((a,b,...),(c,...),z)."""
    Ns = [int(-a) for a in top if a <= 0 and a == int(a)]
    N = min(Ns)                                   # series stops at k=N
    tot = S.Zero
    for k in range(N + 1):
        num = prod([_poch(a, k) for a in top])
        den = prod([_poch(b, k) for b in bot])
        tot += num / den * z**k / factorial(k)
    return tot

beta = symbols('beta', real=True)
fact = factorial

# interior test angles (avoid 0, pi where cos/sin(beta/2) could vanish)
PTS = [Rational(p, 17) * pi for p in (1, 3, 5, 8, 11, 14)]
TOL = S(10) ** (-25)


def jrange(J):
    """M values J, J-1, ..., -J."""
    return [J - k for k in range(int(2 * J) + 1)]


def krange(*constraints):
    """Integer k >= 0 with every constraint(k) >= 0; return the valid list.
    Each constraint is a lambda giving a factorial argument."""
    ks = []
    k = 0
    # upper bound: stop once we've passed all plausible k (2J+2 is safe)
    while k <= 60:
        if all(c(k) >= 0 for c in constraints):
            ks.append(k)
        k += 1
    return ks


def pref(J, M, Mp):
    return sqrt(fact(J + M) * fact(J - M) * fact(J + Mp) * fact(J - Mp))


# ---- eq 4.3.2 ----
def d_432(J, M, Mp):
    ks = krange(lambda k: J - M - k, lambda k: J - Mp - k, lambda k: M + Mp + k)
    s = sum((-1)**k * cos(beta/2)**(M+Mp+2*k) * sin(beta/2)**(2*J-M-Mp-2*k) /
            (fact(k)*fact(J-M-k)*fact(J-Mp-k)*fact(M+Mp+k)) for k in ks)
    return (-1)**(J-Mp) * pref(J, M, Mp) * s


# ---- eq 4.3.3 ----
def d_433(J, M, Mp):
    ks = krange(lambda k: J + M - k, lambda k: J + Mp - k, lambda k: k - M - Mp)
    s = sum((-1)**k * cos(beta/2)**(2*k-M-Mp) * sin(beta/2)**(2*J+M+Mp-2*k) /
            (fact(k)*fact(J+M-k)*fact(J+Mp-k)*fact(k-M-Mp)) for k in ks)
    return (-1)**(J+M) * pref(J, M, Mp) * s      # 4.3.3 phase corrected M'->M


# ---- eq 4.3.4 ----
def d_434(J, M, Mp):
    ks = krange(lambda k: J + M - k, lambda k: J - Mp - k, lambda k: Mp - M + k)
    s = sum((-1)**k * cos(beta/2)**(2*J-2*k+M-Mp) * sin(beta/2)**(2*k-M+Mp) /
            (fact(k)*fact(J+M-k)*fact(J-Mp-k)*fact(Mp-M+k)) for k in ks)
    return pref(J, M, Mp) * s


# ---- eq 4.3.5 ----
def d_435(J, M, Mp):
    ks = krange(lambda k: J - M - k, lambda k: J + Mp - k, lambda k: M - Mp + k)
    s = sum((-1)**k * cos(beta/2)**(2*J-2*k-M+Mp) * sin(beta/2)**(2*k+M-Mp) /
            (fact(k)*fact(J-M-k)*fact(J+Mp-k)*fact(M-Mp+k)) for k in ks)
    return (-1)**(M-Mp) * pref(J, M, Mp) * s


# ---- eq 4.3.6  (general CG form; J1-J2 = M', |J1-J2|<=J<=J1+J2) ----
def d_436(J, M, Mp):
    from sympy.physics.wigner import clebsch_gordan as CG
    J1 = (J + Mp) / 2
    J2 = (J - Mp) / 2                       # J1 - J2 = Mp, J1 + J2 = J
    pre = sqrt(fact(J1 + J2 + J + 1) * fact(J1 + J2 - J) / (2*J + 1))
    s = S.Zero
    for m1 in jrange(J1):
        m2 = M - m1
        if m2 < -J2 or m2 > J2:
            continue
        cg = CG(J1, J2, J, m1, m2, M)   # sympy sig: (j1,j2,j3,m1,m2,m3)
        if cg == 0:
            continue
        s += ((-1)**(J2 + m2) * cg *
              cos(beta/2)**(J1+J2+m1-m2) * sin(beta/2)**(J1+J2-m1+m2) /
              sqrt(fact(J1+m1)*fact(J1-m1)*fact(J2+m2)*fact(J2-m2)))
    return pre * s


TRIG = [("4.3.2", d_432), ("4.3.3", d_433), ("4.3.4", d_434),
        ("4.3.5", d_435), ("4.3.6", d_436)]

x = symbols('x')   # stands for cos(beta) in the Rodrigues forms


def _rodrigues(phase, sqrtarg, plo, phi_, order, alo, ahi, J):
    """Generic differential form:
      phase / 2^J * sqrt(sqrtarg) * (1-x)^plo (1+x)^phi_
        * d^order/dx^order [ (1-x)^alo (1+x)^ahi ]   with x=cos beta."""
    inner = (1 - x)**alo * (1 + x)**ahi
    dexpr = diff(inner, x, int(order))
    full = phase / S(2)**J * sqrt(sqrtarg) * (1 - x)**plo * (1 + x)**phi_ * dexpr
    return full.subs(x, cos(beta))


# ---- eq 4.3.7 .. 4.3.10 (exactly as printed in Chap4.tex) ----
def d_437(J, M, Mp):
    return _rodrigues((-1)**(J-Mp), fact(J+M)/(fact(J-M)*fact(J+Mp)*fact(J-Mp)),
                      (Mp-M)/2, -(M+Mp)/2, J-M, J-Mp, J+Mp, J)

def d_438(J, M, Mp):
    return _rodrigues((-1)**(J+M), fact(J-M)/(fact(J+M)*fact(J+Mp)*fact(J-Mp)),
                      (M-Mp)/2, (M+Mp)/2, J+M, J+Mp, J-Mp, J)

def d_439(J, M, Mp):
    return _rodrigues((-1)**(J-Mp), fact(J+Mp)/(fact(J+M)*fact(J-M)*fact(J-Mp)),
                      (M-Mp)/2, -(M+Mp)/2, J-Mp, J-M, J+M, J)   # inner J-M,J+M

def d_4310(J, M, Mp):
    return _rodrigues((-1)**(J+M), fact(J-Mp)/(fact(J+M)*fact(J-M)*fact(J+Mp)),
                      (Mp-M)/2, (M+Mp)/2, J+Mp, J+M, J-M, J)    # denom (J+M)!, inner J+M,J-M

DIFF = [("4.3.7", d_437), ("4.3.8", d_438), ("4.3.9", d_439), ("4.3.10", d_4310)]


# ---- 4.3.14 / 4.3.15 : mu, nu, s, xi ----
def _mnsx(J, M, Mp):
    mu = abs(M - Mp)
    nu = abs(M + Mp)
    s = J - (mu + nu) / 2
    xi = S(1) if Mp >= M else (-1)**(Mp - M)
    return mu, nu, s, xi


# ---- eq 4.3.13 : Jacobi form ----
def d_4313(J, M, Mp):
    mu, nu, s, xi = _mnsx(J, M, Mp)
    return (xi * sqrt(fact(s)*fact(s+mu+nu)/(fact(s+mu)*fact(s+nu)))
            * sin(beta/2)**mu * cos(beta/2)**nu * jacobi(s, mu, nu, cos(beta)))


# ---- eq 4.3.16 .. 4.3.23 : hypergeometric forms (as printed) ----
def d_4316(J, M, Mp):
    mu, nu, s, xi = _mnsx(J, M, Mp)
    return (xi/fact(mu)*sqrt(fact(s+mu+nu)*fact(s+mu)/(fact(s)*fact(s+nu)))
            * sin(beta/2)**mu * cos(beta/2)**nu
            * hyper((-s, s+mu+nu+1), (mu+1,), sin(beta/2)**2))

def d_4317(J, M, Mp):
    mu, nu, s, xi = _mnsx(J, M, Mp)
    return (xi/fact(mu)*sqrt(fact(s+mu+nu)*fact(s+mu)/(fact(s)*fact(s+nu)))
            * sin(beta/2)**mu * cos(beta/2)**(-nu)
            * hyper((s+mu+1, -s-nu), (mu+1,), sin(beta/2)**2))

def d_4318(J, M, Mp):
    mu, nu, s, xi = _mnsx(J, M, Mp)
    return ((-1)**s*xi/fact(nu)*sqrt(fact(s+mu+nu)*fact(s+nu)/(fact(s)*fact(s+mu)))
            * sin(beta/2)**mu * cos(beta/2)**nu
            * hyper((-s, s+mu+nu+1), (nu+1,), cos(beta/2)**2))

def d_4319(J, M, Mp):
    mu, nu, s, xi = _mnsx(J, M, Mp)
    return ((-1)**s*xi/fact(nu)*sqrt(fact(s+mu+nu)*fact(s+nu)/(fact(s)*fact(s+mu)))
            * sin(beta/2)**(-mu) * cos(beta/2)**nu
            * hyper((s+nu+1, -s-mu), (nu+1,), cos(beta/2)**2))

def d_4320(J, M, Mp):
    mu, nu, s, xi = _mnsx(J, M, Mp)
    return (xi/fact(mu)*sqrt(fact(s+mu+nu)*fact(s+mu)/(fact(s)*fact(s+nu)))
            * sin(beta/2)**mu * cos(beta/2)**(2*s+nu)
            * hyper((-s, -s-nu), (mu+1,), -sin(beta/2)**2/cos(beta/2)**2))

def d_4321(J, M, Mp):
    mu, nu, s, xi = _mnsx(J, M, Mp)
    return ((-1)**s*xi/fact(nu)*sqrt(fact(s+mu+nu)*fact(s+nu)/(fact(s)*fact(s+mu)))
            * sin(beta/2)**(2*s+mu) * cos(beta/2)**nu
            * hyper((-s, -s-mu), (nu+1,), -cos(beta/2)**2/sin(beta/2)**2))

def d_4322(J, M, Mp):
    mu, nu, s, xi = _mnsx(J, M, Mp)
    return ((-1)**s*xi*fact(2*s+mu+nu)
            / sqrt(fact(s)*fact(s+mu+nu)*fact(s+mu)*fact(s+nu))
            * sin(beta/2)**(2*s+mu) * cos(beta/2)**nu
            * hyper((-s, -s-mu), (-2*s-mu-nu,), 1/sin(beta/2)**2))

def d_4323(J, M, Mp):
    mu, nu, s, xi = _mnsx(J, M, Mp)
    return (xi*fact(2*s+mu+nu)
            / sqrt(fact(s)*fact(s+mu+nu)*fact(s+mu)*fact(s+nu))
            * sin(beta/2)**mu * cos(beta/2)**(2*s+nu)
            # source misprint corrected: book has -1/cos^2, correct is +1/cos^2
            * hyper((-s, -s-nu), (-2*s-mu-nu,), 1/cos(beta/2)**2))

HYP = [("4.3.16", d_4316), ("4.3.17", d_4317), ("4.3.18", d_4318),
       ("4.3.19", d_4319), ("4.3.20", d_4320), ("4.3.21", d_4321),
       ("4.3.22", d_4322), ("4.3.23", d_4323)]


def check_pointwise(tag, fn, Jset):
    worst = S.Zero
    nbad = 0
    for J in Jset:
        for M in jrange(J):
            for Mp in jrange(J):
                ref = wigner_d(J, M, Mp)
                got = fn(J, M, Mp)
                for p in PTS:
                    d = Abs((ref - got).subs(beta, p)).evalf(30)
                    if d > worst:
                        worst = d
                    if d > TOL:
                        nbad += 1
                        if nbad <= 4:
                            print(f"    MISMATCH {tag}: J={J} M={M} M'={Mp} "
                                  f"beta={p}  |diff|={float(d):.2e}")
    status = "PASS" if worst < TOL else "FAIL"
    print(f"  eq {tag:8s}  {status}   worst |diff| = {float(worst):.2e}")
    return worst < TOL


def main():
    Jset = [Rational(1, 2), S(1), Rational(3, 2), S(2), Rational(5, 2)]
    print(f"Section 4.3 explicit d-forms vs wigner_d  "
          f"(J up to {Jset[-1]}, {len(PTS)} angles)\n")
    ok = True
    print("Sec 4.3.1  trigonometric sums")
    for tag, fn in TRIG:
        ok &= check_pointwise(tag, fn, Jset)

    # differential + special-function forms: keep J modest (derivatives / series)
    Jsmall = [Rational(1, 2), S(1), Rational(3, 2), S(2)]
    print("\nSec 4.3.2  differential (Rodrigues) representations")
    for tag, fn in DIFF:
        ok &= check_pointwise(tag, fn, Jsmall)
    print("\nSec 4.3.4  Jacobi-polynomial form")
    ok &= check_pointwise("4.3.13", d_4313, Jset)
    print("\nSec 4.3.5  hypergeometric forms")
    for tag, fn in HYP:
        ok &= check_pointwise(tag, fn, Jsmall)

    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
