#!/usr/bin/env python3
r"""
Numerical check of the symmetry relations listed in Section 8.4 of
Varshalovich, Moskalev & Khersonskii, "Quantum Theory of Angular Momentum".

Each relation is tested on a small set of randomly chosen *valid* angular
momenta (a, b, c) and projections (alpha, beta), i.e. triangle inequality
|a-b| <= c <= a+b, integer perimeter a+b+c, and projections in the allowed
ranges.  Configurations are rejected unless the reference 3jm symbol (and
hence the reference Clebsch-Gordan coefficient) is non-zero, so every test
is non-vacuous.

Conventions (matching the book and sympy):
    3jm symbol      ( a b c ; alpha beta gamma )   with  alpha+beta+gamma = 0
                    == sympy.physics.wigner.wigner_3j(a,b,c,alpha,beta,gamma)
    Clebsch-Gordan  C_{a alpha, b beta}^{c gamma} = <a alpha, b beta | c gamma>
                    with gamma = alpha+beta
                    == sympy.physics.wigner.clebsch_gordan(a,b,c,alpha,beta,gamma)

Relations covered
    8.4.2  3jm symbols          eq. (8.4.58) .. (8.4.62)
    8.4.3  Clebsch-Gordan       eq. (8.4.63) .. (8.4.68)
    8.4.5  inversion / t-reversal eq. (8.4.76),(8.4.77)  [reduce to (64)/identity]

Not covered (OCR-damaged in the source and/or requiring negative-j analytic
continuation that sympy does not implement):
    8.4.1  R-symbol relations   eq. (8.4.54)..(8.4.57)
    8.4.4  "mirror" symmetry    eq. (8.4.69)..(8.4.74)

Usage:
    python3 symmetries_8_4.py [--n N] [--seed S] [--jmax J]
"""
from __future__ import annotations

import argparse
import random

from sympy import Rational, sqrt, S
from sympy.physics.wigner import wigner_3j, clebsch_gordan

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
TOL = S(10) ** (-20)          # exact algebraic values -> compare to 20 digits


def _eq(x, y) -> bool:
    """True if the two (exact) sympy expressions are numerically equal."""
    return abs((x - y).evalf(30)) < TOL


def ph(exponent) -> int:
    """(-1)**exponent for an integer-valued exponent."""
    return int(S.NegativeOne ** exponent)


def three_j(a, b, c, al, be, ga):
    return wigner_3j(a, b, c, al, be, ga)


def cg(a, b, c, al, be, ga):
    return clebsch_gordan(a, b, c, al, be, ga)


# ---------------------------------------------------------------------------
# random valid configurations
# ---------------------------------------------------------------------------
def _rand_j(jmax):
    """Random momentum in {0, 1/2, 1, ..., jmax}."""
    return Rational(random.randint(0, int(2 * jmax)), 2)


def random_configs(n, jmax, rng_reject=10000):
    """Return n distinct (a,b,c,alpha,beta) with a non-zero reference 3jm/CG."""
    out, seen, tries = [], set(), 0
    while len(out) < n and tries < rng_reject:
        tries += 1
        a, b = _rand_j(jmax), _rand_j(jmax)
        clo, chi = abs(a - b), a + b
        c = clo + random.randint(0, int(chi - clo))            # step 1, valid
        al = -a + random.randint(0, int(2 * a))                # in {-a,..,a}
        be = -b + random.randint(0, int(2 * b))                # in {-b,..,b}
        ga = al + be                                           # CG projection
        if abs(ga) > c:
            continue
        key = (a, b, c, al, be)
        if key in seen:
            continue
        # non-vacuous: reference coefficient must be non-zero
        if three_j(a, b, c, al, be, -ga) == 0:
            continue
        seen.add(key)
        out.append(key)
    if len(out) < n:
        raise RuntimeError(f"only found {len(out)}/{n} configs; raise --jmax")
    return out


# ---------------------------------------------------------------------------
# the relations
#
# Each entry maps a label -> function(a,b,c,al,be) returning (lhs, rhs).
# For 3jm relations the projection triple is (al, be, ga) with ga = -(al+be).
# For CG  relations the coupling is C_{a al, b be}^{c gam} with gam = al+be.
# ---------------------------------------------------------------------------
def three_j_relations():
    rels = {}

    # eq (8.4.58) -- permutations of columns  (last phase corrected a+b+c;
    # the source prints (-1)^{alpha+b+c}, an obvious OCR typo)
    def f58a(a, b, c, al, be):
        g = -(al + be)
        return three_j(a, b, c, al, be, g), three_j(b, c, a, be, g, al)

    def f58b(a, b, c, al, be):
        g = -(al + be)
        return three_j(a, b, c, al, be, g), three_j(c, a, b, g, al, be)

    def f58c(a, b, c, al, be):
        g = -(al + be)
        return three_j(a, b, c, al, be, g), ph(a + b + c) * three_j(a, c, b, al, g, be)

    def f58d(a, b, c, al, be):
        g = -(al + be)
        return three_j(a, b, c, al, be, g), ph(a + b + c) * three_j(b, a, c, be, al, g)

    def f58e(a, b, c, al, be):
        g = -(al + be)
        return three_j(a, b, c, al, be, g), ph(a + b + c) * three_j(c, b, a, g, be, al)

    # eq (8.4.59) -- change of sign of all projections
    def f59(a, b, c, al, be):
        g = -(al + be)
        return three_j(a, b, c, al, be, g), ph(a + b + c) * three_j(a, b, c, -al, -be, -g)

    # eq (8.4.60) -- Regge, replacement of arguments (first<->third R rows)
    def f60(a, b, c, al, be):
        g = -(al + be)
        H = Rational(1, 2)
        rhs = ph(a + b + c) * three_j(
            (b + c + al) * H, (a + c + be) * H, (a + b + g) * H,
            a - (b + c - al) * H, b - (a + c - be) * H, c - (a + b - g) * H,
        )
        return three_j(a, b, c, al, be, g), rhs

    # eq (8.4.61) -- Regge, transposition of the R-symbol
    def f61(a, b, c, al, be):
        g = -(al + be)
        H = Rational(1, 2)
        rhs = three_j(
            a, (b + c - al) * H, (b + c + al) * H,
            -b + c, (b - c - al) * H - g, (b - c + al) * H + g,
        )
        return three_j(a, b, c, al, be, g), rhs

    # eq (8.4.62) -- the six group representatives (all equal the LHS)
    def _f62(idx):
        def f(a, b, c, al, be):
            g = -(al + be)
            H = Rational(1, 2)
            forms = [
                # form 1
                (a, (b + c - al) * H, (b + c + al) * H,
                 c - b, (b + be - c - g) * H, (b - be - c + g) * H),
                # form 2
                (b, (a + c - be) * H, (a + c + be) * H,
                 a - c, (-a - al + c + g) * H, (-a + al + c - g) * H),
                # form 3
                (c, (a + b - g) * H, (a + b + g) * H,
                 b - a, (a + al - b - be) * H, (a - al - b + be) * H),
                # form 4
                ((b + c - al) * H, (a + c - be) * H, (a + b - g) * H,
                 a - al - (b + c - al) * H, b - be - (a + c - be) * H, c - g - (a + b - g) * H),
                # form 5
                ((b + c + al) * H, (a + c + be) * H, (a + b + g) * H,
                 (b + c + al) * H - a - al, (a + c + be) * H - b - be, (a + b + g) * H - c - g),
            ]
            return three_j(a, b, c, al, be, g), three_j(*forms[idx])
        return f

    rels["eq 8.4.58a  (abc)=(bca)"] = f58a
    rels["eq 8.4.58b  (abc)=(cab)"] = f58b
    rels["eq 8.4.58c  (abc)=+/-(acb)"] = f58c
    rels["eq 8.4.58d  (abc)=+/-(bac)"] = f58d
    rels["eq 8.4.58e  (abc)=+/-(cba)"] = f58e
    rels["eq 8.4.59   projection sign flip"] = f59
    rels["eq 8.4.60   Regge (arg replacement)"] = f60
    rels["eq 8.4.61   Regge (transposition)"] = f61
    for i in range(5):
        rels[f"eq 8.4.62   group repr. #{i + 1}"] = _f62(i)
    return rels


def cg_relations():
    rels = {}

    # eq (8.4.63) -- five CG symmetry relations
    def f63a(a, b, c, al, be):
        g = al + be
        return cg(a, b, c, al, be, g), ph(a + b - c) * cg(b, a, c, be, al, g)

    def f63b(a, b, c, al, be):
        g = al + be
        return cg(a, b, c, al, be, g), ph(a - al) * sqrt((2 * c + 1) / (2 * b + 1)) * cg(a, c, b, al, -g, -be)

    def f63c(a, b, c, al, be):
        g = al + be
        return cg(a, b, c, al, be, g), ph(a - al) * sqrt((2 * c + 1) / (2 * b + 1)) * cg(c, a, b, g, -al, be)

    def f63d(a, b, c, al, be):
        g = al + be
        return cg(a, b, c, al, be, g), ph(b + be) * sqrt((2 * c + 1) / (2 * a + 1)) * cg(c, b, a, -g, be, -al)

    def f63e(a, b, c, al, be):
        g = al + be
        return cg(a, b, c, al, be, g), ph(b + be) * sqrt((2 * c + 1) / (2 * a + 1)) * cg(b, c, a, -be, g, al)

    # eq (8.4.64) -- sign flip of all projections
    def f64(a, b, c, al, be):
        g = al + be
        return cg(a, b, c, al, be, g), ph(a + b - c) * cg(a, b, c, -al, -be, -g)

    # eq (8.4.65)-(8.4.66) -- Regge shuffle, no weight factor
    def f65(a, b, c, al, be):
        g = al + be
        H = Rational(1, 2)
        ap = (a + al + b + be) * H
        alp = (a + al - b - be) * H
        bp = (a - al + b - be) * H
        bep = (a - al - b + be) * H
        cp = c
        gp = a - b
        return cg(a, b, c, al, be, g), cg(ap, bp, cp, alp, bep, gp)

    # eq (8.4.67)-(8.4.68) -- Regge shuffle with weight factor
    def f67(a, b, c, al, be):
        g = al + be
        H = Rational(1, 2)
        ap = (b - be + c + g) * H
        alp = (2 * (a + al) - (b - be) - (c + g)) * H
        bp = (a - al + c + g) * H
        bep = (2 * (b + be) - (a - al) - (c + g)) * H
        cp = (a - al + b - be) * H
        gp = ((a - al) + (b - be) - 2 * (c - g)) * H
        return cg(a, b, c, al, be, g), ph(b + be) * sqrt((2 * c + 1) / (2 * cp + 1)) * cg(ap, bp, cp, alp, bep, gp)

    rels["eq 8.4.63a  swap a<->b"] = f63a
    rels["eq 8.4.63b  a alpha, c -gamma"] = f63b
    rels["eq 8.4.63c  c gamma, a -alpha"] = f63c
    rels["eq 8.4.63d  c -gamma, b beta"] = f63d
    rels["eq 8.4.63e  b -beta, c gamma"] = f63e
    rels["eq 8.4.64   projection sign flip"] = f64
    rels["eq 8.4.65/66 Regge shuffle"] = f65
    rels["eq 8.4.67/68 Regge shuffle (wtd)"] = f67

    # eq (8.4.77) -- time reversal reduces to (8.4.64); include as an explicit
    # check that C(-al,-be,-gam) = (-1)^{a+b-c} C(al,be,gam)
    def f77(a, b, c, al, be):
        g = al + be
        return cg(a, b, c, -al, -be, -g), ph(a + b - c) * cg(a, b, c, al, be, g)

    rels["eq 8.4.77   time reversal (=64)"] = f77
    return rels


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def run(n, seed, jmax):
    random.seed(seed)
    configs = random_configs(n, jmax)

    print(f"Section 8.4 symmetry check -- seed={seed}, jmax={jmax}, "
          f"{len(configs)} random configurations\n")
    print("Configurations  (a, b, c ; alpha, beta, gamma):")
    for (a, b, c, al, be) in configs:
        print(f"    ({a}, {b}, {c} ; {al}, {be}, {al + be})")
    print()

    all_ok = True
    for title, table in (("3jm symbols  (Sec. 8.4.2)", three_j_relations()),
                         ("Clebsch-Gordan  (Sec. 8.4.3 / 8.4.5)", cg_relations())):
        print(f"=== {title} ===")
        for label, fn in table.items():
            npass, first_bad = 0, None
            for cfg in configs:
                lhs, rhs = fn(*cfg)
                if _eq(lhs, rhs):
                    npass += 1
                elif first_bad is None:
                    first_bad = (cfg, lhs, rhs)
            ok = npass == len(configs)
            all_ok &= ok
            mark = "OK  " if ok else "FAIL"
            print(f"  [{mark}] {label:38s} {npass}/{len(configs)}")
            if first_bad:
                cfg, lhs, rhs = first_bad
                print(f"         counterexample {cfg}: lhs={lhs}  rhs={rhs}")
        print()

    print("Not tested (OCR-damaged / negative-j continuation, unsupported):")
    print("    Sec. 8.4.1  R-symbol relations      eq. (8.4.54)-(8.4.57)")
    print("    Sec. 8.4.4  'mirror' symmetry       eq. (8.4.69)-(8.4.74)")
    print()
    print("ALL RELATIONS HOLD" if all_ok else "SOME RELATIONS FAILED -- see above")
    return all_ok


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Check Section 8.4 symmetry relations.")
    p.add_argument("--n", type=int, default=8, help="number of random configs")
    p.add_argument("--seed", type=int, default=20260805, help="RNG seed")
    p.add_argument("--jmax", type=str, default="5/2", help="max momentum, e.g. 5/2")
    args = p.parse_args()
    jmax = Rational(args.jmax)
    ok = run(args.n, args.seed, jmax)
    raise SystemExit(0 if ok else 1)
