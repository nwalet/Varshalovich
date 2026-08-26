#!/usr/bin/env python3
r"""
Checks for Section 4.4 (Symmetries of d^J_{MM'}(beta) and D^J_{MM'}(a,b,g))
of Chapter 4, Varshalovich-Moskalev-Khersonskii.

d^J_{MM'} is the validated VMK helper in wigner_d.py;
D^J_{MM'}(a,b,g) = e^{-iMa} d^J_{MM'}(b) e^{-iM'g}.

Every relation below is checked as LHS-RHS = 0 numerically at several
(a,b,g) tuples over a range of (J,M,M').

  eq 4.4.1   symmetry relations of d (5 rows)
  eq 4.4.4   periodicity of D (4 relations)
  eq 4.4.5   the (a~,g~) re-expression
  eq 4.4.6   D(a+pi, pi-b, -g)      = (-1)^J        D_{M,-M'}(a,b,g)
  eq 4.4.7   D(a-pi, pi-b, pi-g)    = (-1)^{J+M'}   D_{M,-M'}(a,b,g)
  eq 4.4.8   D(a,b,g-pi)           = (-1)^{M'}     D_{M,M'}(a,b,g)

Usage:  python3 check_4_4.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy import I, exp, pi, Rational, S, Abs, symbols
from wigner_d import wigner_d, beta as _BETA

# numeric sample angle tuples (generic, to avoid accidental coincidences)
TRIP = [(Rational(1,5)*pi, Rational(1,4)*pi, Rational(1,7)*pi),
        (Rational(2,7)*pi, Rational(1,3)*pi, Rational(1,6)*pi),
        (Rational(1,3)*pi, Rational(2,5)*pi, Rational(3,8)*pi)]
BETAS = [Rational(1,5)*pi, Rational(1,3)*pi, Rational(2,5)*pi, Rational(3,7)*pi]
TOL = 1e-24


_dcache = {}
def d(J, M, Mp, bb):
    key = (J, M, Mp)
    if key not in _dcache:
        _dcache[key] = wigner_d(J, M, Mp)      # symbolic d^J_{MM'}(beta), once
    return _dcache[key].subs(_BETA, bb)


def Dv(J, M, Mp, aa, bb, gg):
    return exp(-I*M*aa) * d(J, M, Mp, bb) * exp(-I*Mp*gg)


def ph(x):
    """(-1)**x, robust for half-integer x (returns exp(i*pi*x))."""
    return exp(I*pi*x)


def jrange(J):
    return [J - k for k in range(int(2*J) + 1)]


def _report(tag, worst):
    ok = worst < TOL
    print(f"  {tag:52s} {'PASS' if ok else 'FAIL'}   worst={worst:.1e}")
    return ok


# ---------- eq 4.4.1 : d symmetries (test each equality segment) ----------
def check_d_relations(Jset):
    print("eq 4.4.1  symmetry relations of d")
    # each entry: (label, lhs(J,M,Mp,b), rhs(J,M,Mp,b))
    rels = [
        # row 1  (all at +beta)
        ("d = (-1)^{M-M'} d_{-M,-M'}",
         lambda J,M,Mp,b: d(J,M,Mp,b),
         lambda J,M,Mp,b: ph(M-Mp)*d(J,-M,-Mp,b)),
        ("d = (-1)^{M-M'} d_{M',M}   [tex: d_{M'M'}]",
         lambda J,M,Mp,b: d(J,M,Mp,b),
         lambda J,M,Mp,b: ph(M-Mp)*d(J,Mp,M,b)),
        ("d = d_{-M',-M}             [tex: d_{-M'-M'}]",
         lambda J,M,Mp,b: d(J,M,Mp,b),
         lambda J,M,Mp,b: d(J,-Mp,-M,b)),
        # row 2  (LHS at -beta)
        ("d(-b) = (-1)^{M-M'} d(b)",
         lambda J,M,Mp,b: d(J,M,Mp,-b),
         lambda J,M,Mp,b: ph(M-Mp)*d(J,M,Mp,b)),
        ("d(-b) = d_{M',M}(b)",
         lambda J,M,Mp,b: d(J,M,Mp,-b),
         lambda J,M,Mp,b: d(J,Mp,M,b)),
        # row 3  (LHS at pi-beta)
        ("d(pi-b) = (-1)^{J-M'} d_{-M,M'}(b)",
         lambda J,M,Mp,b: d(J,M,Mp,pi-b),
         lambda J,M,Mp,b: ph(J-Mp)*d(J,-M,Mp,b)),
        ("d(pi-b) = (-1)^{J+M} d_{M,-M'}(b)",
         lambda J,M,Mp,b: d(J,M,Mp,pi-b),
         lambda J,M,Mp,b: ph(J+M)*d(J,M,-Mp,b)),
        # row 4  (periodicity 2*pi*n, test n=1 both signs)
        ("d(b+2pi) = (-1)^{2J} d(b)",
         lambda J,M,Mp,b: d(J,M,Mp,b+2*pi),
         lambda J,M,Mp,b: ph(2*J)*d(J,M,Mp,b)),
        # row 5  (b+(2n+1)pi, n=0, +sign)
        ("d(b+pi) = (-1)^{J-M'} d_{M,-M'}(b)",
         lambda J,M,Mp,b: d(J,M,Mp,b+pi),
         lambda J,M,Mp,b: ph(J-Mp)*d(J,M,-Mp,b)),
        ("d(b-pi) = (-1)^{-J-M'} d_{M,-M'}(b)",
         lambda J,M,Mp,b: d(J,M,Mp,b-pi),
         lambda J,M,Mp,b: ph(-J-Mp)*d(J,M,-Mp,b)),
    ]
    ok = True
    for label, lhs, rhs in rels:
        worst = 0.0
        for J in Jset:
            for M in jrange(J):
                for Mp in jrange(J):
                    for b in BETAS:
                        v = Abs((lhs(J,M,Mp,b)-rhs(J,M,Mp,b))).evalf(30)
                        worst = max(worst, float(v))
        ok &= _report(label, worst)
    return ok


# ---------- eqs 4.4.4 - 4.4.8 : D relations ----------
def check_D_relations(Jset):
    print("\neqs 4.4.4-4.4.8  D-function relations")
    n = 1
    rels = [
        # 4.4.4 periodicity
        ("4.4.4a D(a,b+2pi,g)=(-1)^{2J}D(a,b,g)",
         lambda J,M,Mp,a,b,g: Dv(J,M,Mp,a,b+2*pi,g),
         lambda J,M,Mp,a,b,g: ph(2*J)*Dv(J,M,Mp,a,b,g)),
        ("4.4.4b D(a,b+pi,g)=(-1)^{J-M'}D_{M,-M'}(a,b,-g)",
         lambda J,M,Mp,a,b,g: Dv(J,M,Mp,a,b+pi,g),
         lambda J,M,Mp,a,b,g: ph(J-Mp)*Dv(J,M,-Mp,a,b,-g)),
        ("4.4.4c D(a+pi,b,g)=(-i)^{2M}D(a,b,g)",
         lambda J,M,Mp,a,b,g: Dv(J,M,Mp,a+pi,b,g),
         lambda J,M,Mp,a,b,g: (-I)**(2*M)*Dv(J,M,Mp,a,b,g)),
        ("4.4.4d D(a,b,g+pi)=(-i)^{2M'}D(a,b,g)",
         lambda J,M,Mp,a,b,g: Dv(J,M,Mp,a,b,g+pi),
         lambda J,M,Mp,a,b,g: (-I)**(2*Mp)*Dv(J,M,Mp,a,b,g)),
        # 4.4.5
        ("4.4.5 D(a~,b,g~)=e^{iM(a-a~)}D(a,b,g)e^{iM'(g-g~)}",
         lambda J,M,Mp,a,b,g: Dv(J,M,Mp,a/2,b,g/3),
         lambda J,M,Mp,a,b,g: exp(I*M*(a-a/2))*Dv(J,M,Mp,a,b,g)*exp(I*Mp*(g-g/3))),
        # 4.4.6
        ("4.4.6 D(a+pi,pi-b,-g)=(-1)^J D_{M,-M'}(a,b,g)",
         lambda J,M,Mp,a,b,g: Dv(J,M,Mp,a+pi,pi-b,-g),
         lambda J,M,Mp,a,b,g: ph(J)*Dv(J,M,-Mp,a,b,g)),
        # 4.4.7
        ("4.4.7 D(a-pi,pi-b,pi-g)=(-1)^{J+M'}D_{M,-M'}(a,b,g)",
         lambda J,M,Mp,a,b,g: Dv(J,M,Mp,a-pi,pi-b,pi-g),
         lambda J,M,Mp,a,b,g: ph(J+Mp)*Dv(J,M,-Mp,a,b,g)),
        # 4.4.8
        ("4.4.8 D(a,b,g-pi)=(-1)^{M'}D(a,b,g)",
         lambda J,M,Mp,a,b,g: Dv(J,M,Mp,a,b,g-pi),
         lambda J,M,Mp,a,b,g: ph(Mp)*Dv(J,M,Mp,a,b,g)),
    ]
    ok = True
    for label, lhs, rhs in rels:
        worst = 0.0
        for J in Jset:
            for M in jrange(J):
                for Mp in jrange(J):
                    for (a,b,g) in TRIP:
                        v = Abs((lhs(J,M,Mp,a,b,g)-rhs(J,M,Mp,a,b,g))).evalf(30)
                        worst = max(worst, float(v))
        ok &= _report(label, worst)
    return ok


# ---------- eq 4.4.2 : the big D symmetry array (transcribed from the scan) ----------
def check_D_array(Jset):
    print("\neq 4.4.2  D symmetry array (31 relations, from scan)")
    # index selectors and prefactors
    mm = lambda M, Mp: (M, Mp);   MM = lambda M, Mp: (-M, -Mp)
    pm = lambda M, Mp: (Mp, M);   PM = lambda M, Mp: (-Mp, -M)
    P1 = lambda e, h: 1; Pe = lambda e, h: e; Ph = lambda e, h: h; Peh = lambda e, h: e*h

    def table(a, b, g):
        return [
            # row 1 (base is col-0): angles (a,b,g)
            (Peh, MM, (a, b, g), False), (Ph, mm, (a, b, g), True), (Pe, MM, (a, b, g), True),
            # row 2 (g,b,a)
            (Pe, pm, (g, b, a), False), (Ph, PM, (g, b, a), False), (Peh, pm, (g, b, a), True), (P1, PM, (g, b, a), True),
            # row 3 (a,-b,g)
            (Pe, mm, (a, -b, g), False), (Ph, MM, (a, -b, g), False), (Peh, mm, (a, -b, g), True), (P1, MM, (a, -b, g), True),
            # row 4 (g,-b,a)
            (P1, pm, (g, -b, a), False), (Peh, PM, (g, -b, a), False), (Ph, pm, (g, -b, a), True), (Pe, PM, (g, -b, a), True),
            # row 5 (-a,b,-g)
            (Pe, MM, (-a, b, -g), False), (Ph, mm, (-a, b, -g), False), (Peh, MM, (-a, b, -g), True), (P1, mm, (-a, b, -g), True),
            # row 6 (-g,b,-a)
            (P1, PM, (-g, b, -a), False), (Peh, pm, (-g, b, -a), False), (Ph, PM, (-g, b, -a), True), (Pe, pm, (-g, b, -a), True),
            # row 7 (-a,-b,-g)
            (P1, MM, (-a, -b, -g), False), (Peh, mm, (-a, -b, -g), False), (Ph, MM, (-a, -b, -g), True), (Pe, mm, (-a, -b, -g), True),
            # row 8 (-g,-b,-a)
            (Pe, PM, (-g, -b, -a), False), (Ph, pm, (-g, -b, -a), False), (Peh, PM, (-g, -b, -a), True), (P1, pm, (-g, -b, -a), True),
        ]
    worst = 0.0
    worst_idx = None
    for J in Jset:
        for M in jrange(J):
            for Mp in jrange(J):
                for (a, b, g) in TRIP[:2]:
                    base = complex(Dv(J, M, Mp, a, b, g).evalf(30))
                    e = complex(ph(Mp - M).evalf(30))
                    h = complex((exp(-I*2*M*a - I*2*Mp*g)).evalf(30))
                    for i, (pref, idx, ang, conj) in enumerate(table(a, b, g)):
                        Mi, Mip = idx(M, Mp)
                        val = complex(Dv(J, Mi, Mip, *ang).evalf(30))
                        if conj:
                            val = val.conjugate()
                        dev = abs(base - pref(e, h) * val)
                        if dev > worst:
                            worst, worst_idx = dev, (i, J, M, Mp)
    ok = worst < 1e-10   # double-precision (complex()) floor is ~1e-16
    print(f"  all 31 relations                                    "
          f"{'PASS' if ok else 'FAIL'}   worst={worst:.1e}"
          + ("" if ok else f"  at entry {worst_idx}"))
    return ok


def main():
    Jset = [Rational(1,2), S(1), Rational(3,2), S(2)]
    print(f"Section 4.4 symmetries  (J up to {Jset[-1]})\n")
    ok = check_d_relations(Jset)
    ok &= check_D_array(Jset)
    ok &= check_D_relations(Jset)
    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
