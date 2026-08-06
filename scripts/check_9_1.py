#!/usr/bin/env python3
r"""
Checks for Section 9.1 (Definition) of Chapter 9, Varshalovich, Moskalev &
Khersonskii -- 6j symbols, Racah coefficients and the R-symbol.

Numeric (sympy angular-momentum functions):
  eq 9.1.8   sum of four CG = 6j                       (definition)
  eq 9.1.9   6j orthogonality (sum over j12) = delta_{j23 j23'}
  eq 9.1.10  6j orthogonality (sum over j23) = delta_{j12 j12'}   [OCR: RHS delta index]
  eq 9.1.11  6j = (-1)^{a+b+d+e} W(a b e d; c f)       (Racah phase)

Symbolic (R-symbol algebra, a..f free):
  eq 9.1.13  R_{i alpha} in terms of a..f
  eq 9.1.14  inverse relations 2a=R13+R24=R14+R23, ...
  eq 9.1.16  row/column sums of R
  eq 9.1.17/18  parametrisation R_{i alpha}=A_i - B_alpha
  eq 9.1.19  sum A_i = sum B_alpha = 2(a+..+f)
  eq 9.1.20  inverse relations via A_i, B_alpha

Usage:  python3 check_9_1.py
"""
from sympy import Rational, S, sqrt, symbols, simplify
from sympy.physics.wigner import clebsch_gordan as CG, wigner_6j, racah

H = Rational(1, 2)


def tri(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) == int(a + b + c)


def proj(j):
    return [-j + i for i in range(int(2 * j) + 1)]


def C(a, b, c, al, be, ga):
    if a < 0 or b < 0 or c < 0 or abs(al) > a or abs(be) > b or abs(ga) > c:
        return S.Zero
    if not tri(a, b, c) or al + be != ga:
        return S.Zero
    return CG(a, b, c, al, be, ga)


def w6(a, b, c, d, e, f):
    """6j, 0 when any of its four triads violates the triangle rule."""
    if not (tri(a, b, c) and tri(a, e, f) and tri(d, b, f) and tri(d, e, c)):
        return S.Zero
    return wigner_6j(a, b, c, d, e, f)


# --------------------------------------------------------------------------
# eq 9.1.8 : sum of four Clebsch-Gordan coefficients = 6j
# --------------------------------------------------------------------------
def eq_8(j1, j2, j3, j12, j23, j, m):
    """LHS sum over m1,m2,m3,m12,m23 (m'=m, j'=j) vs RHS."""
    tot = S.Zero
    for m1 in proj(j1):
        for m2 in proj(j2):
            for m3 in proj(j3):
                m12 = m1 + m2
                m23 = m2 + m3
                tot += (C(j12, j3, j, m12, m3, m) * C(j1, j2, j12, m1, m2, m12)
                        * C(j1, j23, j, m1, m23, m) * C(j2, j3, j23, m2, m3, m23))
    rhs = (S.NegativeOne ** (j1 + j2 + j3 + j)
           * sqrt((2 * j12 + 1) * (2 * j23 + 1)) * w6(j1, j2, j12, j3, j, j23))
    return simplify(tot - rhs) == 0


# --------------------------------------------------------------------------
# eq 9.1.9 / 9.1.10 : orthogonality
# --------------------------------------------------------------------------
def eq_9(j1, j2, j3, j, j23, j23p, jmax):
    tot = S.Zero
    for j12 in [Rational(i, 2) for i in range(0, int(2 * jmax) + 1)]:
        tot += ((2 * j12 + 1) * (2 * j23 + 1)
                * w6(j1, j2, j12, j3, j, j23) * w6(j1, j2, j12, j3, j, j23p))
    return simplify(tot - (1 if j23 == j23p else 0)) == 0


def eq_10(j1, j2, j3, j, j12, j12p, jmax):
    tot = S.Zero
    for j23 in [Rational(i, 2) for i in range(0, int(2 * jmax) + 1)]:
        tot += ((2 * j12 + 1) * (2 * j23 + 1)
                * w6(j1, j2, j12, j3, j, j23) * w6(j1, j2, j12p, j3, j, j23))
    return simplify(tot - (1 if j12 == j12p else 0)) == 0


# --------------------------------------------------------------------------
# eq 9.1.11 : 6j = (-1)^{a+b+d+e} W(a b e d; c f)
#   W(a b e d; c f) = racah(a, b, e, d, c, f) in sympy's argument order
# --------------------------------------------------------------------------
def eq_11(a, b, c, d, e, f):
    lhs = w6(a, b, c, d, e, f)
    rhs = S.NegativeOne ** (a + b + d + e) * racah(a, b, e, d, c, f)
    return simplify(lhs - rhs) == 0


# --------------------------------------------------------------------------
# symbolic R-symbol algebra (eq 9.1.13 - 9.1.20)
# --------------------------------------------------------------------------
def r_symbol_algebra():
    a, b, c, d, e, f = symbols('a b c d e f')
    R = {
        (1, 1): -c + d + e, (1, 2): b + d - f, (1, 3): a + e - f, (1, 4): a + b - c,
        (2, 1): -b + d + f, (2, 2): c + d - e, (2, 3): a - b + c, (2, 4): a - e + f,
        (3, 1): -a + e + f, (3, 2): -a + b + c, (3, 3): c - d + e, (3, 4): b - d + f,
    }
    out = []

    # eq 9.1.14 : inverse relations
    inv = [
        (2 * a, R[1, 3] + R[2, 4]), (2 * a, R[1, 4] + R[2, 3]),
        (2 * d, R[1, 1] + R[2, 2]), (2 * d, R[1, 2] + R[2, 1]),
        (2 * b, R[1, 2] + R[3, 4]), (2 * b, R[1, 4] + R[3, 2]),
        (2 * e, R[1, 1] + R[3, 3]), (2 * e, R[1, 3] + R[3, 1]),
        (2 * c, R[2, 2] + R[3, 3]), (2 * c, R[2, 3] + R[3, 2]),
        (2 * f, R[2, 1] + R[3, 4]), (2 * f, R[2, 4] + R[3, 1]),
    ]
    out.append(("eq 9.1.14  inverse relations",
                all(simplify(l - r) == 0 for l, r in inv)))

    # eq 9.1.16 : column sums and grand total
    colsum = [
        (sum(R[i, 1] for i in (1, 2, 3)), 2 * (d + e + f) - a - b - c),
        (sum(R[i, 3] for i in (1, 2, 3)), 2 * (a + c + e) - b - d - f),
        (sum(R[i, 2] for i in (1, 2, 3)), 2 * (b + c + d) - a - e - f),
        (sum(R[i, 4] for i in (1, 2, 3)), 2 * (a + b + f) - c - d - e),
        (sum(R.values()), 2 * (a + b + c + d + e + f)),
    ]
    out.append(("eq 9.1.16  row/column sums",
                all(simplify(l - r) == 0 for l, r in colsum)))

    # eq 9.1.17/18 : R_{i alpha} = A_i - B_alpha
    A = {1: a + b + d + e, 2: a + c + d + f, 3: b + c + e + f}
    B = {1: a + b + c, 2: a + e + f, 3: b + d + f, 4: c + d + e}
    out.append(("eq 9.1.17/18  R_{i,al}=A_i-B_al",
                all(simplify(A[i] - B[al] - R[i, al]) == 0 for i, al in R)))

    # eq 9.1.19 : sum A_i = sum B_alpha = 2(a+..+f)
    tot = 2 * (a + b + c + d + e + f)
    out.append(("eq 9.1.19  sum A_i = sum B_al = 2(a+..+f)",
                simplify(sum(A.values()) - tot) == 0
                and simplify(sum(B.values()) - tot) == 0))

    # eq 9.1.20 : inverse relations via A_i, B_alpha
    inv2 = [
        (2 * a, A[1] + A[2] - B[3] - B[4]), (2 * d, A[1] + A[2] - B[1] - B[2]),
        (2 * b, A[1] + A[3] - B[2] - B[4]), (2 * e, A[1] + A[3] - B[1] - B[3]),
        (2 * c, A[2] + A[3] - B[2] - B[3]), (2 * f, A[2] + A[3] - B[1] - B[4]),
    ]
    out.append(("eq 9.1.20  inverse relations via A,B",
                all(simplify(l - r) == 0 for l, r in inv2)))
    return out


def run():
    print("Section 9.1 checks\n")
    ok = True

    # numeric cases (small, triangle-valid)
    # (j1, j2, j3, j12, j23, j, m) with all four recoupling triads valid
    cases8 = [(1, 1, 1, 1, 1, 1, 0),
              (1, 1, 1, 2, 2, 1, 0),
              (H, 1, H, H, H, 1, 0),
              (Rational(3, 2), 1, H, H, Rational(3, 2), 1, 0)]
    r8 = all(eq_8(*x) for x in cases8)
    print(f"  [{'OK  ' if r8 else 'FAIL'}] eq 9.1.8   four-CG sum = 6j        ({len(cases8)} cases)")
    ok &= r8

    r9 = eq_9(1, 1, 1, 1, 1, 1, 3) and eq_9(1, 1, 1, 1, 1, 2, 3) and eq_9(Rational(3,2), 1, H, 1, H, Rational(3,2), 3)
    print(f"  [{'OK  ' if r9 else 'FAIL'}] eq 9.1.9   orthogonality (sum j12)")
    ok &= r9

    r10 = eq_10(1, 1, 1, 1, 1, 1, 3) and eq_10(1, 1, 1, 1, 1, 2, 3) and eq_10(Rational(3,2), 1, H, 1, H, Rational(3,2), 3)
    print(f"  [{'OK  ' if r10 else 'FAIL'}] eq 9.1.10  orthogonality (sum j23)  [confirms RHS = delta_{{j12,j12'}}]")
    ok &= r10

    cases11 = [(1, 1, 1, 1, 1, 1), (Rational(3,2), 1, H, 1, Rational(3,2), 1),
               (2, 1, 1, 1, 2, 1), (H, H, 1, H, H, 1)]
    r11 = all(eq_11(*x) for x in cases11)
    print(f"  [{'OK  ' if r11 else 'FAIL'}] eq 9.1.11  6j = (-1)^(a+b+d+e) W    ({len(cases11)} cases)")
    ok &= r11

    print()
    for label, res in r_symbol_algebra():
        print(f"  [{'OK  ' if res else 'FAIL'}] {label}")
        ok &= res

    print("\nALL 9.1 CHECKS PASS" if ok else "\nSOME 9.1 CHECKS FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
