#!/usr/bin/env python3
r"""
Checks for Section 9.2 of Chapter 9 (Varshalovich, Moskalev & Khersonskii):
general expressions for the 6j symbols.

Numeric (compared against sympy's wigner_6j / wigner_3j):
  eq 9.2.1, 9.2.2   Racah single-sum formulas
  eq 9.2.3 - 9.2.6  other single-sum formulas
  eq 9.2.8 - 9.2.12 hypergeometric 4F3(1) representations
  eq 9.2.13         6j as a sum of four 3jm symbols (phase check)

Not checked here: eq 9.2.7 (Bargmann, R-symbol double sum) and the
quasi-binomial forms eq 9.2.14-9.2.24 (structural).

Usage:  python3 check_9_2.py
"""
from sympy import Rational, S, sqrt, factorial as fac, RisingFactorial as rf, simplify
from sympy.physics.wigner import wigner_6j, wigner_3j

H = Rational(1, 2)


def tri(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) == int(a + b + c)


def valid6(a, b, c, d, e, f):
    # the four triads of {a b c / d e f}
    return tri(a, b, c) and tri(c, d, e) and tri(a, e, f) and tri(b, d, f)


def w6(a, b, c, d, e, f):
    return wigner_6j(a, b, c, d, e, f) if valid6(a, b, c, d, e, f) else S.Zero


def w3(a, b, c, al, be, ga):
    if abs(al) > a or abs(be) > b or abs(ga) > c or al + be + ga != 0:
        return S.Zero
    if not tri(a, b, c):
        return S.Zero
    return wigner_3j(a, b, c, al, be, ga)


def D(a, b, c):
    return sqrt(fac(a + b - c) * fac(a - b + c) * fac(-a + b + c) / fac(a + b + c + 1))


def proj(j):
    return [-j + i for i in range(int(2 * j) + 1)]


def term(sign, num_args, den_args):
    """sign * prod(num!)/prod(den!), or 0 if any factorial argument is negative."""
    if any(x < 0 for x in num_args + den_args):
        return S.Zero
    num = S.One
    for x in num_args:
        num *= fac(x)
    den = S.One
    for x in den_args:
        den *= fac(x)
    return sign * num / den


def nmax(a, b, c, d, e, f):
    return int(a + b + c + d + e + f) + 3


# --------------------------------------------------------------------------
# single-sum formulas
# --------------------------------------------------------------------------
def eq1(a, b, c, d, e, f):
    # CORRECTED: (a+c+d+f-n)! and (b+c+e+f-n)! belong in the DENOMINATOR
    # (the tex prints them in the numerator -- OCR error).
    s = sum(term((-1) ** n, [n + 1],
                 [n - a - b - c, n - c - d - e, n - a - e - f, n - b - d - f,
                  a + b + d + e - n, a + c + d + f - n, b + c + e + f - n])
            for n in range(nmax(a, b, c, d, e, f) + 1))
    return D(a, b, c) * D(c, d, e) * D(a, e, f) * D(b, d, f) * s


def eq2(a, b, c, d, e, f):
    # CORRECTED: (-a+c-d+f+n)! and (-b+c-e+f+n)! belong in the DENOMINATOR.
    s = sum(term((-1) ** n, [a + b + d + e + 1 - n],
                 [n, a + b - c - n, -c + d + e - n, a + e - f - n, b + d - f - n,
                  -a + c - d + f + n, -b + c - e + f + n])
            for n in range(nmax(a, b, c, d, e, f) + 1))
    return (-1) ** (a + b + d + e) * D(a, b, c) * D(c, d, e) * D(a, e, f) * D(b, d, f) * s


def eq3(a, b, c, d, e, f):
    s = sum(term((-1) ** n, [-a + b + c + n, c - d + e + n, a - c + d + f - n],
                 [n, a - e + f - n, -b + d + f - n, -a + b - d + e + n, b + c + e - f + 1 + n])
            for n in range(nmax(a, b, c, d, e, f) + 1))
    return (-1) ** (a + c + d + f) * D(a, e, f) * D(b, d, f) / (D(a, b, c) * D(c, d, e)) * s


def eq4(a, b, c, d, e, f):
    s = sum(term((-1) ** n, [a - b + d + e - n, -b + c + e + f - n, a + c + d + f + 1 - n],
                 [n, a - b + c - n, -b + d + f - n, a + e + f + 1 - n, c + d + e + 1 - n])
            for n in range(nmax(a, b, c, d, e, f) + 1))
    return (-1) ** (a + c + d + f) * D(a, b, c) * D(b, d, f) / (D(a, e, f) * D(c, d, e)) * s


def eq5(a, b, c, d, e, f):
    s = sum(term((-1) ** n, [-a + e + f + n, b - d + f + n, a + c + d - f - n],
                 [n, a + e - f - n, b + d - f - n, -a + c - d + f + n, 2 * f + 1 + n])
            for n in range(nmax(a, b, c, d, e, f) + 1))
    pref = ((-1) ** (a + b + d + e) * D(a, b, c) * D(c, d, e) * D(a, e, f) * D(b, d, f)
            * fac(a + e + f + 1) * fac(b + d + f + 1)
            / (fac(a + b - c) * fac(a - b + c) * fac(-c + d + e) * fac(c + d - e)
               * fac(-a + e + f) * fac(b - d + f)))
    return pref * s


def eq6(a, b, c, d, e, f):
    s = sum(term((-1) ** n, [2 * b - n, b + c - e + f - n, b + c + e + f + 1 - n],
                 [n, -a + b + c - n, b - d + f - n, a + b + c + 1 - n, b + d + f + 1 - n])
            for n in range(nmax(a, b, c, d, e, f) + 1))
    pref = ((-1) ** (b + c + e + f) * D(a, b, c) * D(c, d, e) * D(a, e, f) * D(b, d, f)
            * fac(a + b + c + 1) * fac(b + d + f + 1)
            / (fac(a + b - c) * fac(c - d + e) * fac(c + d - e) * fac(a - e + f)
               * fac(-a + e + f) * fac(b + d - f)))
    return pref * s


# --------------------------------------------------------------------------
# hypergeometric 4F3(1) forms
# --------------------------------------------------------------------------
def F43(top, bot):
    negs = [int(-t) for t in top if t == int(t) and t <= 0]
    N = min(negs)
    tot = S.Zero
    for k in range(0, N + 1):
        num = S.One
        for t in top:
            num *= rf(t, k)
        den = fac(k)
        for b in bot:
            den *= rf(b, k)
        tot += num / den
    return tot


def eq8(a, b, c, d, e, f):
    pref = ((-1) ** (a + b + d + e) * D(a, b, c) * D(c, d, e) * D(a, e, f) * D(b, d, f)
            * fac(a + b + d + e + 1)
            / (fac(a + b - c) * fac(-c + d + e) * fac(a + e - f) * fac(b + d - f)
               * fac(-a + c - d + f) * fac(-b + c - e + f)))
    top = [-a - b + c, c - d - e, -a - e + f, -b - d + f]
    bot = [-a - b - d - e - 1, -a + c - d + f + 1, -b + c - e + f + 1]
    return pref * F43(top, bot)


def eq9(a, b, c, d, e, f):
    pref = ((-1) ** (a + c + d + f) * D(a, e, f) * D(b, d, f) / (D(a, b, c) * D(c, d, e))
            * fac(-a + b + c) * fac(c - d + e) * fac(a - c + d + f)
            / (fac(a - e + f) * fac(-b + d + f) * fac(-a + b - d + e) * fac(b + c + e - f + 1)))
    top = [-a + b + c + 1, c - d + e + 1, -a + e - f, b - d - f]
    bot = [-a + c - d - f, -a + b - d + e + 1, b + c + e - f + 2]
    return pref * F43(top, bot)


def eq10(a, b, c, d, e, f):
    pref = ((-1) ** (a + c + d + f) * D(a, b, c) * D(b, d, f) / (D(a, e, f) * D(c, d, e))
            * fac(a - b + d + e) * fac(-b + c + e + f) * fac(a + c + d + f + 1)
            / (fac(a - b + c) * fac(-b + d + f) * fac(a + e + f + 1) * fac(c + d + e + 1)))
    top = [-a + b - c, b - d - f, -a - e - f - 1, -c - d - e - 1]
    bot = [-a + b - d - e, b - c - e - f, -a - c - d - f - 1]
    return pref * F43(top, bot)


def eq11(a, b, c, d, e, f):
    pref = ((-1) ** (a + b + d + e) * D(a, b, c) * D(c, d, e) * D(a, e, f) * D(b, d, f)
            * fac(a + e + f + 1) * fac(b + d + f + 1) * fac(a + c + d - f)
            / (fac(a + b - c) * fac(a - b + c) * fac(-c + d + e) * fac(c + d - e)
               * fac(a + e - f) * fac(b + d - f) * fac(-a + c - d + f) * fac(2 * f + 1)))
    top = [-a - e + f, -b - d + f, -a + e + f + 1, b - d + f + 1]
    bot = [-a - c - d + f, -a + c - d + f + 1, 2 * f + 2]
    return pref * F43(top, bot)


def eq12(a, b, c, d, e, f):
    pref = ((-1) ** (b + c + e + f) * D(a, b, c) * D(c, d, e) * D(a, e, f) * D(b, d, f)
            * fac(2 * b) * fac(b + c - e + f) * fac(b + c + e + f + 1)
            / (fac(-a + b + c) * fac(a + b - c) * fac(c - d + e) * fac(c + d - e)
               * fac(a - e + f) * fac(-a + e + f) * fac(b + d - f) * fac(b - d + f)))
    top = [a - b - c, -b + d - f, -a - b - c - 1, -b - d - f - 1]
    bot = [-2 * b, -b - c + e - f, -b - c - e - f - 1]
    return pref * F43(top, bot)


# --------------------------------------------------------------------------
# eq 9.2.13 : 6j as sum of four 3jm  (phase check)
# --------------------------------------------------------------------------
def eq13(a, b, c, d, e, f, use_eps=True):
    tot = S.Zero
    for al in proj(a):
        for be in proj(b):
            ga = -al - be
            if abs(ga) > c:
                continue
            for ep in proj(e):
                ph = al + ep          # phi  (from a,e,f: al+ep-phi=0)
                if abs(ph) > f:
                    continue
                de = be + ph          # delta (from d,b,f: -de+be+ph=0)
                if abs(de) > d:
                    continue
                # exponent: literal text is d+e+f+delta+e+phi; ε-reading uses epsilon
                expo = d + e + f + de + (ep if use_eps else e) + ph
                tot += ((-1) ** expo * w3(a, b, c, al, be, ga) * w3(a, e, f, al, ep, -ph)
                        * w3(d, b, f, -de, be, ph) * w3(d, e, c, de, -ep, ga))
    return tot


CASES = [(1, 1, 1, 1, 1, 1), (2, 1, 1, 1, 2, 1), (Rational(3, 2), 1, H, 1, Rational(3, 2), 1),
         (2, 2, 2, 1, 1, 2), (H, H, 1, H, H, 1), (2, Rational(3, 2), H, 1, Rational(3, 2), 1),
         (2, 2, 2, 2, 2, 2), (2, 2, 1, 2, 2, 1), (3, 2, 1, 2, 2, 2), (2, 2, 2, 1, 2, 1)]


import math

TOL = 1e-15


def _num(v):
    """Python float value, or None if symbolic / non-finite (nan/inf)."""
    if v is None:
        return None
    try:
        d = complex(v.evalf(30))
    except Exception:
        return None
    if not (math.isfinite(d.real) and math.isfinite(d.imag)) or abs(d.imag) > 1e-12:
        return None
    return d.real


def check(name, fn):
    # some closed forms have degenerate points (a prefactor factorial of a
    # negative argument, or a 4F3 bottom parameter = nonpositive integer);
    # they represent the 6j only where finite, so skip non-finite evaluations.
    ok = True
    tested = 0
    for x in CASES:
        if not valid6(*x):
            continue
        v = _num(fn(*x))
        if v is None:
            continue
        tested += 1
        if abs(v - _num(w6(*x))) > TOL:
            ok = False
            break
    ok = ok and tested > 0
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:28s} ({tested} well-defined cases)")
    return ok


def run():
    print("Section 9.2 checks\n")
    ok = True
    for name, fn in [("eq 9.2.1  Racah sum", eq1), ("eq 9.2.2  Racah sum (n->..)", eq2),
                     ("eq 9.2.3  single sum", eq3), ("eq 9.2.4  single sum", eq4),
                     ("eq 9.2.5  single sum", eq5), ("eq 9.2.6  single sum", eq6),
                     ("eq 9.2.8  4F3(1)", eq8), ("eq 9.2.9  4F3(1)", eq9),
                     ("eq 9.2.10 4F3(1)", eq10), ("eq 9.2.11 4F3(1)", eq11),
                     ("eq 9.2.12 4F3(1)", eq12)]:
        ok &= check(name, fn)

    # eq 9.2.13 : compare the two phase readings
    eps_ok = all(abs(_num(eq13(*x, use_eps=True)) - _num(w6(*x))) <= TOL for x in CASES if valid6(*x))
    lit_ok = all(abs(_num(eq13(*x, use_eps=False)) - _num(w6(*x))) <= TOL for x in CASES if valid6(*x))
    print(f"  [{'OK  ' if eps_ok else 'FAIL'}] eq 9.2.13 four-3jm sum, phase (-1)^(d+e+f+delta+epsilon+phi)")
    print(f"  [{'--  ' if not lit_ok else '??  '}] eq 9.2.13 literal phase (...+e+...) holds: {lit_ok}")
    ok &= eps_ok

    print("\nALL 9.2 CHECKS PASS" if ok else "\nSOME 9.2 CHECKS FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
