#!/usr/bin/env python3
r"""
Checks for Section 9.4 (symmetries) of Chapter 9, Varshalovich, Moskalev &
Khersonskii.

  eq 9.4.2  classical symmetries of the 6j symbol (24 forms)
  eq 9.4.3  Regge symmetries of the 6j symbol (with s1,s2,s3)
  eq 9.4.5  classical symmetries of the Racah W-coefficient (24 forms)
  eq 9.4.7  Regge symmetries of the Racah W-coefficient

Each listed form is evaluated for several valid argument sets and compared with
the reference value; a mismatch flags an OCR error in the printed arguments.
(eq 9.4.1 R-symbol row/column permutation and eq 9.4.8/9.4.10 "mirror"
j -> -j-1 continuation are structural, not covered here.)

Usage:  python3 check_9_4.py
"""
from sympy import Rational, S
from sympy.physics.wigner import wigner_6j, racah

H = Rational(1, 2)


def tri(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) == int(a + b + c)


def valid6(a, b, c, d, e, f):
    return tri(a, b, c) and tri(c, d, e) and tri(a, e, f) and tri(b, d, f)


def w6(a, b, c, d, e, f):
    return wigner_6j(a, b, c, d, e, f) if valid6(a, b, c, d, e, f) else S.Zero


# VMK W(a b e d ; c f) with {a b c / d e f} = (-1)^{a+b+d+e} W(a b e d; c f);
# in sympy the string W(p q r s; t u) is racah(p, q, r, s, t, u).
def Wv(p, q, r, s, t, u):
    # racah needs valid triads; guard by mapping back to the 6j it represents
    # W(p q r s; t u) ~ 6j {p q t / s r u}; validity:
    if not (tri(p, q, t) and tri(s, r, t) and tri(p, r, u) and tri(s, q, u)):
        return S.Zero
    return racah(p, q, r, s, t, u)


# reference 6j argument order
def ref6(a, b, c, d, e, f):
    return w6(a, b, c, d, e, f)


# ---- eq 9.4.2 : 24 classical 6j forms (each equals {a b c / d e f}) ----
def forms_942(a, b, c, d, e, f):
    return [
        (a, b, c, d, e, f), (a, c, b, d, f, e), (b, a, c, e, d, f), (b, c, a, e, f, d),
        (c, a, b, f, d, e), (c, b, a, f, e, d),
        (a, e, f, d, b, c), (a, f, e, d, c, b), (e, a, f, b, d, c), (e, f, a, b, c, d),
        (f, a, e, c, d, b), (f, e, a, c, b, d),
        (d, e, c, a, b, f), (d, c, e, a, f, b), (e, d, c, b, a, f), (e, c, d, b, f, a),
        (c, d, e, f, a, b), (c, e, d, f, b, a),
        (d, b, f, a, e, c), (d, f, b, a, c, e), (b, d, f, e, a, c), (b, f, d, e, c, a),
        (f, d, b, c, a, e), (f, b, d, c, e, a),
    ]


# ---- eq 9.4.3 : Regge 6j symmetries ----
def forms_943(a, b, c, d, e, f):
    s1 = Rational(b + c + e + f, 2)
    s2 = Rational(a + c + d + f, 2)
    s3 = Rational(a + b + d + e, 2)
    return [
        (a, b, c, d, e, f),
        (a, s1 - b, s1 - c, d, s1 - e, s1 - f),
        (s2 - a, b, s2 - c, s2 - d, e, s2 - f),
        (s3 - a, s3 - b, c, s3 - d, s3 - e, f),
        (s2 - d, s3 - e, s1 - f, s2 - a, s3 - b, s1 - c),
        (s3 - d, s1 - e, s2 - f, s3 - a, s1 - b, s2 - c),
    ]


# ---- eq 9.4.5 : 24 classical Racah forms; (phase_key, args) ----
def forms_945(a, b, c, d, e, f):
    e1 = (-1) ** (b + e - c - f)
    e2 = (-1) ** (a + d - c - f)
    L = [
        (1, (a, b, e, d, c, f)), (1, (d, e, b, a, c, f)), (1, (e, d, a, b, c, f)), (1, (b, a, d, e, c, f)),
        (1, (a, e, b, d, f, c)), (1, (d, b, e, a, f, c)), (1, (b, d, a, e, f, c)), (1, (e, a, d, b, f, c)),
        (e1, (a, c, f, d, b, e)), (e1, (d, f, c, a, b, e)), (e1, (f, d, a, c, b, e)), (e1, (c, a, d, f, b, e)),
        (e1, (a, f, c, d, e, b)), (e1, (d, c, f, a, e, b)), (e1, (c, d, a, f, e, b)), (e1, (f, a, d, c, e, b)),
        (e2, (c, b, e, f, a, d)), (e2, (f, e, b, c, a, d)), (e2, (e, f, c, b, a, d)), (e2, (b, c, f, e, a, d)),
        (e2, (c, e, b, f, d, a)), (e2, (f, b, e, c, d, a)), (e2, (b, f, c, e, d, a)), (e2, (e, c, f, b, d, a)),
    ]
    return L


# ---- eq 9.4.7 : Regge Racah symmetries ----
def forms_947(a, b, c, d, e, f):
    s1 = Rational(b + c + e + f, 2)
    s2 = Rational(a + c + d + f, 2)
    s3 = Rational(a + b + d + e, 2)
    e1 = (-1) ** (b + e - c - f)
    e2 = (-1) ** (a + d - c - f)
    return [
        (1, (a, b, e, d, c, f)),
        (1, (s3 - a, s3 - b, s3 - e, s3 - d, c, f)),
        (e1, (a, s1 - b, s1 - e, d, s1 - c, s1 - f)),
        (e1, (s2 - d, s3 - e, s3 - b, s2 - a, s1 - f, s1 - c)),
        (e2, (s2 - a, b, e, s2 - d, s2 - c, s2 - f)),
        (e2, (s3 - d, s1 - e, s1 - b, s3 - a, s2 - f, s2 - c)),
    ]


CASES = [(2, 2, 2, 2, 2, 2), (3, 2, 1, 2, 2, 2), (2, 2, 2, 1, 1, 2),
         (3, 3, 2, 2, 2, 2), (2, Rational(3, 2), H, 1, Rational(3, 2), 1),
         (Rational(5, 2), 2, H, 2, Rational(3, 2), 1)]


def check6(name, formfn):
    ok = True
    n = 0
    for x in CASES:
        if not valid6(*x):
            continue
        r = ref6(*x)
        for g in formfn(*x):
            n += 1
            if w6(*g) != r:
                ok = False
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:40s} ({n} comparisons)")
    return ok


def _eq(u, v):
    return abs(complex((u - v).evalf(30))) < 1e-18


def checkW(name, formfn):
    ok = True
    n = 0
    for x in CASES:
        if not valid6(*x):
            continue
        ref = Wv(x[0], x[1], x[4], x[3], x[2], x[5])   # W(a b e d; c f)
        for ph, args in formfn(*x):
            n += 1
            if not _eq(S(ph) * Wv(*args), ref):   # S(ph): keep phase exact (-1**-1 -> float in py)
                ok = False
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:40s} ({n} comparisons)")
    return ok


def run():
    print("Section 9.4 symmetry checks\n")
    ok = True
    ok &= check6("eq 9.4.2  classical 6j symmetries", forms_942)
    ok &= check6("eq 9.4.3  Regge 6j symmetries", forms_943)
    ok &= checkW("eq 9.4.5  classical Racah symmetries", forms_945)
    ok &= checkW("eq 9.4.7  Regge Racah symmetries", forms_947)
    print("\nALL 9.4 CHECKS PASS" if ok else "\nSOME 9.4 CHECKS FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
