#!/usr/bin/env python3
r"""
Checks for Section 9.5 (explicit 6j values for special arguments) of Chapter 9,
Varshalovich, Moskalev & Khersonskii.

Covers the closed forms 9.5.1, 9.5.3-9.5.20 (compared to sympy wigner_6j) and
the V_c(a,f,b) system 9.5.21-9.5.38 (V_c reconstructed from the 6j and compared
to the explicit polynomials / special values).

Two source oddities are examined:
  * eq 9.5.21 phase (-1)^{a+b+c+J}: J is undefined -- we test J = f.
  * eq 9.5.7 \sixj{a}{b}{a+b}{a}{e}{-1}: the last argument "-1" looks corrupt;
    we test the hypothesis f = a+b-1.

Usage:  python3 check_9_5.py
"""
import math
from sympy import Rational, S, sqrt, factorial as fac
from sympy.physics.wigner import wigner_6j, racah

H = Rational(1, 2)


def tri(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) == int(a + b + c)


def valid6(a, b, c, d, e, f):
    return tri(a, b, c) and tri(c, d, e) and tri(a, e, f) and tri(b, d, f)


def w6(a, b, c, d, e, f):
    return wigner_6j(a, b, c, d, e, f) if valid6(a, b, c, d, e, f) else S.Zero


def num(v):
    if v is None:
        return None
    try:
        d = complex(v.evalf(30))
    except Exception:
        return None
    if not (math.isfinite(d.real) and math.isfinite(d.imag)) or abs(d.imag) > 1e-12:
        return None
    return d.real


def F(x):
    return fac(x)


def rt(numf, denf):
    """sqrt(prod numf! / prod denf!) or None if any argument negative."""
    if any(x < 0 for x in numf + denf):
        return None
    n = S.One
    for x in numf:
        n *= F(x)
    dd = S.One
    for x in denf:
        dd *= F(x)
    return sqrt(n / dd)


def close(u, v):
    a, b = num(u), num(v)
    return a is not None and b is not None and abs(a - b) < 1e-15


# ----- eq 9.5.1 : one argument zero (6 forms of the 6j) -----
def check_951():
    # each form: (which slot is 0, 6j-arg builder from free params, phase, sqrt-args)
    ok = True
    n = 0
    vals = [H, 1, Rational(3, 2), 2, Rational(5, 2)]
    for p in vals:
        for q in vals:
            for r in vals:
                # 1: {0 b c;d e f}, b=c, e=f  -> {0 q q; r p p}
                tests = [
                    ((0, q, q, r, p, p), (-1) ** (q + p + r) / sqrt((2 * q + 1) * (2 * p + 1))),      # (-1)^{b+e+d}
                    # 2: {a b c;0 e f}, b=f, c=e -> {p q r; 0 r q}
                    ((p, q, r, 0, r, q), (-1) ** (p + q + r) / sqrt((2 * q + 1) * (2 * r + 1))),      # (-1)^{a+b+e}, e=r
                    # 3: {a 0 c;d e f}, a=c, d=f -> {p 0 p; r q r}
                    ((p, 0, p, r, q, r), (-1) ** (p + r + q) / sqrt((2 * p + 1) * (2 * r + 1))),      # (-1)^{a+d+e}
                    # 4: {a b c;d 0 f}, a=f, c=d -> {p q r; r 0 p}
                    ((p, q, r, r, 0, p), (-1) ** (p + q + r) / sqrt((2 * p + 1) * (2 * r + 1))),      # (-1)^{a+b+d}, d=r
                    # 5: {a b 0;d e f}, a=b, d=e -> {p p 0; r r q}
                    ((p, p, 0, r, r, q), (-1) ** (p + r + q) / sqrt((2 * p + 1) * (2 * r + 1))),      # (-1)^{a+e+f}, e=r,f=q
                    # 6: {a b c;d e 0}, a=e, b=d -> {p q r; q p 0}
                    ((p, q, r, q, p, 0), (-1) ** (p + q + r) / sqrt((2 * p + 1) * (2 * q + 1))),      # (-1)^{a+b+c}
                ]
                for args, rhs in tests:
                    if not valid6(*args):
                        continue
                    n += 1
                    if not close(w6(*args), rhs):
                        ok = False
    print(f"  [{'OK  ' if ok else 'FAIL'}] eq 9.5.1  one argument zero        ({n} cases)")
    return ok


# ----- helper: iterate valid parameter tuples for the {a b a+b; d e f} family
def params_sum(kmax=3):
    vals = [H, 1, Rational(3, 2), 2, Rational(5, 2)]
    for a in vals:
        for b in vals:
            for d in vals:
                for e in vals:
                    for f in vals:
                        yield a, b, d, e, f


# ----- eq 9.5.3 : {a b a+b; d e f} -----
def eq3(a, b, d, e, f):
    r = rt([2 * a, 2 * b, a + b + d + e + 1, a + b - d + e, a + b + d - e, -a + e + f, -b + d + f],
           [2 * a + 2 * b + 1, -a - b + d + e, a + e - f, a - e + f, a + e + f + 1, b + d - f, b - d + f, b + d + f + 1])
    return None if r is None else (-1) ** (a + b + d + e) * r


# ----- eq 9.5.8 : {a b a+b; a b f} -----
def eq8(a, b, f):
    if any(x < 0 for x in [a + b - f, a + b + f + 1]):
        return None
    return (-1) ** (2 * a + 2 * b) * F(2 * a) * F(2 * b) / (F(a + b - f) * F(a + b + f + 1))


# ----- eq 9.5.9 : {a b a+b; b a f} -----
def eq9(a, b, f):
    r = rt([], [2 * a - f, 2 * a + f + 1, 2 * b - f, 2 * b + f + 1])
    return None if r is None else (-1) ** (2 * a + 2 * b) * F(2 * a) * F(2 * b) * r


# ----- eq 9.5.10 : {a b a+b; a b a+b} -----
def eq10(a, b):
    return (-1) ** (2 * a + 2 * b) * F(2 * a) * F(2 * b) / F(2 * a + 2 * b + 1)


# ----- eq 9.5.15 : {a b a+b-1; d e f} -----
def eq15(a, b, d, e, f):
    r = rt([2 * a - 1, 2 * b - 1, a + b + d + e, a + b - d + e - 1, a + b + d - e - 1, -a + e + f, -b + d + f],
           [2 * a + 2 * b, -a - b + d + e + 1, a + e - f, a - e + f, a + e + f + 1, b + d - f, b - d + f, b + d + f + 1])
    if r is None:
        return None
    br = 2 * (a * b * (a + b) + (a + b) * f * (f + 1) - a * d * (d + 1) - b * e * (e + 1))
    return (-1) ** (a + b + d + e) * br * r


# ----- eq 9.5.18 : {a b a+b-1; a b f} -----
def eq18(a, b, f):
    if any(x < 0 for x in [a + b - f, a + b + f + 1]):
        return None
    br = 2 * (a * b * (a + b) + (a + b) * f * (f + 1) - a ** 2 * (a + 1) - b ** 2 * (b + 1))
    return (-1) ** (2 * a + 2 * b) * br * F(2 * a - 1) * F(2 * b - 1) / (F(a + b - f) * F(a + b + f + 1))


# ----- eq 9.5.19 : {a b a+b-1; b a f} -----
def eq19(a, b, f):
    r = rt([], [2 * a - f, 2 * a + f + 1, 2 * b - f, 2 * b + f + 1])
    if r is None:
        return None
    br = 2 * ((a + b) * f * (f + 1) - 2 * a * b)
    return (-1) ** (2 * a + 2 * b) * br * F(2 * a - 1) * F(2 * b - 1) * r


def check_closed():
    ok = True
    counts = {}

    def cmp(tag, val, six):
        counts[tag] = counts.get(tag, 0)
        v = num(val)
        if v is None:
            return
        counts[tag] += 1
        if not close(val, six):
            counts[tag] = -abs(counts[tag]) - 1  # mark failure

    for a, b, d, e, f in params_sum():
        c = a + b
        if valid6(a, b, c, d, e, f):
            cmp("eq 9.5.3 ", eq3(a, b, d, e, f), w6(a, b, c, d, e, f))
        # eq8: {a b a+b; a b f}
        if valid6(a, b, a + b, a, b, f):
            cmp("eq 9.5.8 ", eq8(a, b, f), w6(a, b, a + b, a, b, f))
        if valid6(a, b, a + b, b, a, f):
            cmp("eq 9.5.9 ", eq9(a, b, f), w6(a, b, a + b, b, a, f))
        if valid6(a, b, a + b, a, b, a + b):
            cmp("eq 9.5.10", eq10(a, b), w6(a, b, a + b, a, b, a + b))
        # 9.5.15 family c=a+b-1
        if a + b - 1 >= 0 and valid6(a, b, a + b - 1, d, e, f):
            cmp("eq 9.5.15", eq15(a, b, d, e, f), w6(a, b, a + b - 1, d, e, f))
        if a + b - 1 >= 0 and valid6(a, b, a + b - 1, a, b, f):
            cmp("eq 9.5.18", eq18(a, b, f), w6(a, b, a + b - 1, a, b, f))
        if a + b - 1 >= 0 and valid6(a, b, a + b - 1, b, a, f):
            cmp("eq 9.5.19", eq19(a, b, f), w6(a, b, a + b - 1, b, a, f))

    for tag in ["eq 9.5.3 ", "eq 9.5.8 ", "eq 9.5.9 ", "eq 9.5.10", "eq 9.5.15", "eq 9.5.18", "eq 9.5.19"]:
        c = counts.get(tag, 0)
        good = c >= 1
        ok &= good
        print(f"  [{'OK  ' if good else 'FAIL'}] {tag}  ({c if c>=0 else 'FAILED'} cases)")
    return ok


# ----- eq 9.5.21 + V_c system -----
def Vfrom6j(a, f, b, c, J_is_f=True):
    r = rt([2 * a - c, 2 * b - c], [2 * a + c + 1, 2 * b + c + 1])
    if r is None:
        return None
    Jphase = f if J_is_f else 0
    ph = (-1) ** (a + b + c + Jphase)
    denom = ph * r
    return w6(a, a, c, b, b, f) / denom


def Vexplicit(c, a, f, b):
    ta, tb = a * (a + 1), b * (b + 1)
    x = f * (f + 1) - a * (a + 1) - b * (b + 1)
    if c == 0:
        return S.One
    if c == 1:
        return -2 * x
    if c == 2:
        return 6 * x ** 2 + 6 * x - 8 * ta * tb
    if c == 3:
        return -20 * x ** 3 - 80 * x ** 2 - 16 * x * (3 + ta + tb - 3 * ta * tb) + 80 * ta * tb
    if c == 4:
        return (70 * x ** 4 + 700 * x ** 3 + 40 * x ** 2 * (39 + 5 * ta + 5 * tb - 6 * ta * tb)
                + 80 * x * (9 + 6 * ta + 6 * tb - 17 * ta * tb) - 48 * ta * tb * (27 + 4 * ta + 4 * tb - 2 * ta * tb))
    return None


def check_V():
    ok = True
    # confirm J = f by matching V_0..V_4 polynomials
    for J_is_f in (True, False):
        good = 0
        bad = 0
        for a in [1, Rational(3, 2), 2, Rational(5, 2)]:
            for b in [1, Rational(3, 2), 2]:
                for f in [i * H for i in range(0, 11)]:
                    for c in range(0, 5):
                        if not valid6(a, a, c, b, b, f):
                            continue
                        Vc = Vfrom6j(a, f, b, c, J_is_f)
                        Ve = Vexplicit(c, a, f, b)
                        if Vc is None or Ve is None:
                            continue
                        if close(Vc, Ve):
                            good += 1
                        else:
                            bad += 1
        tag = "J=f" if J_is_f else "J=0"
        print(f"        V_0..V_4 with {tag}: {good} match, {bad} mismatch")
        if J_is_f:
            ok &= (bad == 0 and good > 0)

    # special-f forms 9.5.33-9.5.38 (with J=f)
    def V33(c, a, b): return F(2 * b) * F(2 * a + c + 1) / (F(2 * a + 1) * F(2 * b - c))          # f=a-b, a>=b
    def V34(c, a, b): return F(2 * a) * F(2 * b + c + 1) / (F(2 * b + 1) * F(2 * a - c))          # f=b-a, a<=b
    def V37(c, a, b): return (-1) ** (c + 1) * 2 * ((a + b) * c * (c + 1) - 2 * a * b) * F(2 * a - 1) * F(2 * b - 1) / (F(2 * a - c) * F(2 * b - c))  # f=a+b-1
    def V38(c, a, b): return (-1) ** c * F(2 * a) * F(2 * b) / (F(2 * a - c) * F(2 * b - c))       # f=a+b

    sp = [("9.5.33", lambda c, a, b: (a - b, V33(c, a, b)) if a >= b else None),
          ("9.5.34", lambda c, a, b: (b - a, V34(c, a, b)) if a <= b else None),
          ("9.5.37", lambda c, a, b: (a + b - 1, V37(c, a, b))),
          ("9.5.38", lambda c, a, b: (a + b, V38(c, a, b)))]
    for tag, gen in sp:
        good = bad = 0
        for a in [1, Rational(3, 2), 2, Rational(5, 2)]:
            for b in [1, Rational(3, 2), 2]:
                for c in range(0, int(2 * min(a, b)) + 1):
                    g = gen(c, a, b)
                    if g is None:
                        continue
                    f, Vval = g
                    if f < 0 or not valid6(a, a, c, b, b, f):
                        continue
                    Vc = Vfrom6j(a, f, b, c, True)
                    if Vc is None:
                        continue
                    if close(Vc, Vval):
                        good += 1
                    else:
                        bad += 1
        okk = bad == 0 and good > 0
        ok &= okk
        print(f"  [{'OK  ' if okk else 'FAIL'}] eq {tag}  V_c special f      ({good} ok, {bad} bad)")
    return ok


# ----- eq 9.5.7 hypothesis: f = a+b-1 -----
def check_957():
    # RHS as printed
    def rhs(a, b, e):
        r = rt([], [2 * a + 2 * b + 1, -b + e + 1]) if False else None
        inner = rt([2 * a, 2 * b, b + e, 2 * a + b - e, 2 * a + b + e + 1] if False else [], [])
        # (2a-1)!(b+e-1)!/((2a+2b)!(-b+e)!) * sqrt(2a*2b*(b+e)*(2a+b-e)*(2a+b+e+1)/((2a+2b+1)*(-b+e+1)))
        if any(x < 0 for x in [2 * a - 1, b + e - 1, 2 * a + 2 * b, -b + e]):
            return None
        pref = (-1) ** (2 * a + b + e) * F(2 * a - 1) * F(b + e - 1) / (F(2 * a + 2 * b) * F(-b + e))
        rad = Rational(2 * a * 2 * b * (b + e) * (2 * a + b - e) * (2 * a + b + e + 1), 1) / ((2 * a + 2 * b + 1) * (-b + e + 1))
        if rad < 0:
            return None
        return pref * sqrt(rad)
    good = bad = 0
    for a in [1, Rational(3, 2), 2]:
        for b in [1, Rational(3, 2), 2]:
            for e in [i * H for i in range(1, 9)]:
                r = rhs(a, b, e)
                if r is None:
                    continue
                six = w6(a, b, a + b, a, e, a + b - 1)     # hypothesis f=a+b-1
                if not valid6(a, b, a + b, a, e, a + b - 1):
                    continue
                if close(r, six):
                    good += 1
                else:
                    bad += 1
    print(f"  [{'OK  ' if good > 0 and bad == 0 else 'FAIL'}] eq 9.5.7  hypothesis f=a+b-1  ({good} ok, {bad} bad)")
    return good > 0 and bad == 0


def run():
    print("Section 9.5 checks\n")
    ok = True
    ok &= check_951()
    ok &= check_closed()
    print()
    ok &= check_V()
    print()
    ok &= check_957()
    print("\nALL 9.5 CHECKS PASS" if ok else "\nSOME 9.5 CHECKS FAILED / need inspection")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
